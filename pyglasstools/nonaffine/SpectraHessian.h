#ifndef __SPECTRAHESSIAN_H__
#define __SPECTRAHESSIAN_H__

#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>
#include <memory>
#include <chrono>

#include <omp.h>

#include <pyglasstools/ParticleSystem.h>
#include <pyglasstools/potential/PairPotential.h>
#include <pyglasstools/MPIFile.h>
#include "NonAffineManager.h"

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <Aboria.h>
namespace abr = Aboria;

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/StdVector>

typedef Eigen::Triplet<double> Tripletd;
typedef Eigen::SparseMatrix<double> SparseMatd;

#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

class PYBIND11_EXPORT SpectraHessian
{
    protected:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< PairPotential > m_potential;
        std::shared_ptr< HessianManager > m_manager;
        std::shared_ptr< MPI::Communicator > m_comm;
        double max_rcut;
        unsigned int hessian_length;
        
        Eigen::VectorXd eigenvals;
        Eigen::MatrixXd eigenvecs;
        Eigen::MatrixXd nonaffinetensor;
        Eigen::MatrixXd pseudoinverse;
        Eigen::MatrixXd misforce;
        SparseMatd hessian;
        unsigned int nconv; 
    public:

        #pragma omp declare reduction( + : Eigen::MatrixXd : omp_out += omp_in ) initializer( omp_priv = Eigen::MatrixXd::Zero(omp_orig.rows(),omp_orig.cols()) )
        
        SpectraHessian(    std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential,
                    std::shared_ptr< HessianManager > manager, std::shared_ptr< MPI::Communicator > comm);
        virtual ~SpectraHessian(){};
        
        void setSystemData(std::shared_ptr<ParticleSystem> newsysdata)
        {
            m_sysdata = newsysdata;
        }

        void buildHessianandMisForce();
        
        void getEigenPairs();
        void getAllEigenPairs();
        
        std::tuple<bool,double> getEigenvalue(unsigned int index);
        void saveEigenvalue(unsigned int index, std::shared_ptr<MPI::LogFile> logfile);
        
        std::tuple<bool,Eigen::VectorXd> getEigenvector(unsigned int index);
        void saveEigenvector(unsigned int index, std::shared_ptr<MPI::LogFile> logfile);
        
        void calculateNonAffineTensor(); 
        std::tuple<bool,double> getNonAffineTensor(int i, int j); 
        void saveNonAffineTensor(unsigned int i, unsigned int j, std::shared_ptr<MPI::LogFile> logfile);

};

SpectraHessian::SpectraHessian(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential,
                    std::shared_ptr< HessianManager > manager, std::shared_ptr< MPI::Communicator > comm)
            : m_sysdata(sysdata), m_potential(potential), m_manager(manager), m_comm(comm)
{
    if (m_comm->getSizeGlobal() > 1)
    {
        m_manager->notice(5) << "[WARNING] Performing MPI run. The hessian class supports no MPI parallelization on its eigendecomposition." << std::endl;
        m_manager->notice(5) << "[WARNING] MPI parallelization of eigendecomposition only exists in 'slepc' mode." << std::endl;
    }
    double maxforce_rcut = potential->scaled_rcut;
    double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                            std::end(abr::get<diameter>(m_sysdata->particles)) );
    maxforce_rcut *= maxdiameter;
    max_rcut = std::max(maxforce_rcut,maxdiameter); //whichever is maximum
    
    //Construct Hessian
    hessian_length = (unsigned int)m_sysdata->simbox->dim*m_sysdata->particles.size();
    
    //Set size for eigenvectors and eigenvalues
    eigenvecs.resize(hessian_length,1);
    eigenvecs.setZero();
    //Set size for eigenvectors and eigenvalues
    eigenvals.resize(1);
    eigenvals.setZero();
    //Construct pseudoinverse
    pseudoinverse = Eigen::MatrixXd::Zero(hessian_length,hessian_length);
    buildHessianandMisForce();  
    nonaffinetensor = Eigen::MatrixXd::Zero(3,3);
    nconv = 0;
}

void SpectraHessian::buildHessianandMisForce()
{
    unsigned int estimate_nonzero_entries = (unsigned int)m_sysdata->simbox->dim*hessian_length*6;
    std::vector< Tripletd > hessian_triplet;
    hessian_triplet.reserve(estimate_nonzero_entries);
    hessian.resize(hessian_length,hessian_length);
    
    //Construct mismatch force vectors
    misforce.resize(hessian_length, (unsigned int)m_sysdata->simbox->dim*((unsigned int)m_sysdata->simbox->dim+1)/2);
    misforce.setZero();
    for( auto p_i = m_sysdata->particles.begin(); p_i != m_sysdata->particles.end(); ++p_i)
    {
        
        for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                    abr::get<position>(*p_i), max_rcut); p_j != false; ++p_j)
        {
            //Compute a list of local obsercavles 
            //Next we loop through the j-th particles for the virial stress
            //set the distance between particle i and particle j
            Eigen::Vector3d rij(p_j.dx()[0], p_j.dx()[1], p_j.dx()[2]);
            Eigen::Vector3d nij = rij.normalized();
            m_potential->rij = rij;
            
            //Don't forget to set diameters of the potential
            m_potential->di =  abr::get<diameter>(*p_i);
            m_potential->dj =  abr::get<diameter>(*p_j);
            
            int id_i = abr::get<abr::id>(*p_i);
            int id_j = abr::get<abr::id>(*p_j);

            //Make sure the particle is unique
            if (id_i != id_j && m_potential->getRcut() > rij.dot(rij))
            {
                double factor = m_potential->getBondStiffness()+m_potential->getPairForce();
                Eigen::Matrix3d offdiag_ij = -factor*nij*nij.transpose()
                                             +Eigen::Matrix3d::Identity()*m_potential->getPairForce();

                misforce(2*id_i,0) += -factor*rij[0]*nij[0]*nij[0]; 
                misforce(2*id_i+1,0) += -factor*rij[0]*nij[0]*nij[1]; 
                
                misforce(2*id_i,1) += -factor*rij[1]*nij[1]*nij[0]; 
                misforce(2*id_i+1,1) += -factor*rij[1]*nij[1]*nij[1]; 
                
                misforce(2*id_i,2) += -factor*rij[0]*nij[1]*nij[0]; 
                misforce(2*id_i+1,2) += -factor*rij[0]*nij[1]*nij[1]; 
                
                hessian_triplet.push_back(Tripletd( 2*id_i, 2*id_j, offdiag_ij(0,0) ));
                hessian_triplet.push_back(Tripletd( 2*id_i, 2*id_j+1, offdiag_ij(0,1) ));
                hessian_triplet.push_back(Tripletd( 2*id_i+1, 2*id_j, offdiag_ij(1,0) ));
                hessian_triplet.push_back(Tripletd( 2*id_i+1, 2*id_j+1, offdiag_ij(1,1) ));
                
                
                hessian_triplet.push_back(Tripletd( 2*id_i, 2*id_i, -offdiag_ij(0,0) ));
                hessian_triplet.push_back(Tripletd( 2*id_i, 2*id_i+1, -offdiag_ij(0,1) ));
                hessian_triplet.push_back(Tripletd( 2*id_i+1, 2*id_i, -offdiag_ij(1,0) ));
                hessian_triplet.push_back(Tripletd( 2*id_i+1, 2*id_i+1, -offdiag_ij(1,1) ));
            }
        }
    }
    hessian.setFromTriplets(hessian_triplet.begin(), hessian_triplet.end());
}

void SpectraHessian::getEigenPairs()
{
    Eigen::MatrixXd globaleigenvecs;//(hessia
    Eigen::VectorXd globaleigenvals;//(hessia
    double tol = m_manager->lowerbound_tol;
    std::string selrule = m_manager->selrule;
    int nev = m_manager->nev;
    int ncv = m_manager->ncv;
    int maxiter = m_manager->maxiter;

    //only do EigenPairs in the root process
    if (m_comm->isRoot())
    {
        // Construct matrix operation object using the wrapper class SparseGenMatProd
        int info = Spectra::NOT_COMPUTED;
        if (selrule == "LM")
        {
            Spectra::SparseGenMatProd<double> op(hessian);
            Spectra::SymEigsSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double> > eigs(&op, nev, ncv);
            eigs.init();
            nconv = eigs.compute(maxiter, tol,Spectra::SMALLEST_ALGE);
            info = eigs.info();
            globaleigenvals = eigs.eigenvalues().real();
            globaleigenvecs = eigs.eigenvectors().real();
        }    
        else if (selrule == "LR")
        {
            Spectra::SparseGenMatProd<double> op(hessian);
            Spectra::SymEigsSolver< double, Spectra::LARGEST_REAL, Spectra::SparseGenMatProd<double> > eigs(&op, nev, ncv);
            eigs.init();
            nconv = eigs.compute(maxiter, tol,Spectra::SMALLEST_ALGE);
            info = eigs.info();
            globaleigenvals = eigs.eigenvalues().real();
            globaleigenvecs = eigs.eigenvectors().real();
        }    
        else if (selrule == "LA")
        {
            Spectra::SparseGenMatProd<double> op(hessian);
            Spectra::SymEigsSolver< double, Spectra::LARGEST_ALGE, Spectra::SparseGenMatProd<double> > eigs(&op, nev, ncv);
            eigs.init();
            nconv = eigs.compute(maxiter, tol,Spectra::SMALLEST_ALGE);
            info = eigs.info();
            globaleigenvals = eigs.eigenvalues().real();
            globaleigenvecs = eigs.eigenvectors().real();
        }    
        else if (selrule == "SM")
        {
            Spectra::SparseSymShiftSolve<double> op(hessian);
            Spectra::SymEigsShiftSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<double> > eigs(&op, nev, ncv, tol);
            eigs.init();
            nconv = eigs.compute(maxiter, tol,Spectra::SMALLEST_ALGE);
            info = eigs.info();
            globaleigenvals = eigs.eigenvalues().real();
            globaleigenvecs = eigs.eigenvectors().real();
        }    
        else if (selrule == "SR")
        {
            Spectra::SparseSymShiftSolve<double> op(hessian);
            Spectra::SymEigsShiftSolver< double, Spectra::LARGEST_REAL, Spectra::SparseSymShiftSolve<double> > eigs(&op, nev, ncv, tol);
            eigs.init();
            nconv = eigs.compute(maxiter, tol,Spectra::SMALLEST_ALGE);
            info = eigs.info();
            globaleigenvals = eigs.eigenvalues().real();
            globaleigenvecs = eigs.eigenvectors().real();
        }    
        else if (selrule == "SA")
        {
            Spectra::SparseSymShiftSolve<double> op(hessian);
            Spectra::SymEigsShiftSolver< double, Spectra::LARGEST_ALGE, Spectra::SparseSymShiftSolve<double> > eigs(&op, nev, ncv, tol);
            eigs.init();
            nconv = eigs.compute(maxiter, tol,Spectra::SMALLEST_ALGE);
            info = eigs.info();
            globaleigenvals = eigs.eigenvalues().real();
            globaleigenvecs = eigs.eigenvectors().real();
        }    
        else 
        {
            //Do a large MPI ABort
            m_manager->notice(5) << "[ERROR] Selection Rule unrecognized" << std::endl;
            m_comm->abort(1);
        }
    
        if (info == Spectra::SUCCESSFUL)             
        {
            m_manager->notice(5) << "Eigenvalues and Eigenvectors computed successfully." << std::endl;;
            m_manager->notice(5) << "Number of requested eigenpairs:" << nev << std::endl;
            m_manager->notice(5) << "Number of converged eigenpairs:" << nconv << std::endl;
        }
    }
    m_manager->notice(5) << "Distributing results to all child processes" << std::endl;;
    m_comm->bcast(nconv,0); 
    //Once I've figured out scatter with Eigen objects, I'll change this to a single scatter_v command
    div_t divres = div((int)nconv, (int)m_comm->getSizeGlobal());
    std::vector<int> counts(m_comm->getSizeGlobal(), divres.quot);
    for(int i = 0; i < m_comm->getSizeGlobal(); ++i)
    {
        if (i < divres.rem) 
            counts[i] += 1;
        m_manager->notice(5) << counts[i] << std::endl;
    }
    if(m_comm->isRoot())
    {
        eigenvecs = globaleigenvecs.block(0,0,hessian_length,counts[0]);
        eigenvals = globaleigenvals.segment(0,counts[0]);
    }
    for(unsigned int i = 1; i < m_comm->getSizeGlobal(); ++i)
    {
        int starts = std::accumulate(counts.begin(), counts.begin()+i, 0);
        if(m_comm->isRoot())
        {
            Eigen::MatrixXd temp = globaleigenvecs.block(0,starts,hessian_length,counts[i]);
            m_comm->send(temp, i, 0);

            Eigen::VectorXd tempeigval = globaleigenvals.segment(starts,counts[i]);
            m_comm->send(tempeigval, i,1);
        }
        else if(m_comm->getRank() == i)
        {
            m_comm->recv(eigenvecs, 0, 0);
            m_comm->recv(eigenvals, 0, 1);
        }
        else
        {
            continue;
        }
    } 

    m_manager->widenotice(5) << "the size of my eigenvector matrix is: " << eigenvecs.rows() << "x" << eigenvecs.cols() << std::endl;
    m_manager->widenotice(5) << "the size of my eigenvalue vector is: " << eigenvals.rows() << "x" << eigenvals.cols() << std::endl;
    m_comm->barrier();
    //m_manager->widenotice(5) << "the size of my eigenvalue vector is: " << eigenvalue.rows() << "x" << eigenvalue.cols() << std::endl;
}

void SpectraHessian::getAllEigenPairs()
{
    Eigen::MatrixXd globaleigenvecs;//(hessia
    Eigen::VectorXd globaleigenvals;//(hessia
    double tol = m_manager->lowerbound_tol;
    int nev = hessian.rows()-2;
    int ncv = hessian.rows();
    int maxiter = m_manager->maxiter;

    //only do EigenPairs in the root process
    if (m_comm->isRoot())
    {
        // Construct matrix operation object using the wrapper class SparseGenMatProd
        int info = Spectra::NOT_COMPUTED;
        Spectra::SparseGenMatProd<double> op(hessian);
        Spectra::SymEigsSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double> > eigs(&op, nev, ncv);
        eigs.init();
        nconv = eigs.compute(maxiter, tol,Spectra::SMALLEST_ALGE);
        info = eigs.info();
        globaleigenvals = eigs.eigenvalues().real();
        globaleigenvecs = eigs.eigenvectors().real();
    
        if (info == Spectra::SUCCESSFUL)             
        {
            m_manager->notice(5) << "Eigenvalues and Eigenvectors computed successfully." << std::endl;;
            m_manager->notice(5) << "Number of requested eigenpairs:" << nev << std::endl;
            m_manager->notice(5) << "Number of converged eigenpairs:" << nconv << std::endl;
        }
    }
    m_manager->notice(5) << "Distributing results to all child processes" << std::endl;;
    m_comm->bcast(nconv,0); 
     
    //Once I've figured out scatter with Eigen objects, I'll change this to a single scatter_v command
    div_t divres = div((int)nconv, (int)m_comm->getSizeGlobal());
    std::vector<int> counts(m_comm->getSizeGlobal(), divres.quot);
    for(int i = 0; i < m_comm->getSizeGlobal(); ++i)
    {
        if (i < divres.rem) 
            counts[i] += 1;
        m_manager->notice(5) << counts[i] << std::endl;
    }
    if(m_comm->isRoot())
    {
        eigenvecs = globaleigenvecs.block(0,0,hessian_length,counts[0]);
        eigenvals = globaleigenvals.segment(0,counts[0]);
    }
    for(unsigned int i = 1; i < m_comm->getSizeGlobal(); ++i)
    {
        int starts = std::accumulate(counts.begin(), counts.begin()+i, 0);
        if(m_comm->isRoot())
        {
            Eigen::MatrixXd temp = globaleigenvecs.block(0,starts,hessian_length,counts[i]);
            m_comm->send(temp, i, 0);

            Eigen::VectorXd tempeigval = globaleigenvals.segment(starts,counts[i]);
            m_comm->send(tempeigval, i,1);
        }
        else if(m_comm->getRank() == i)
        {
            m_comm->recv(eigenvecs, 0, 0);
            m_comm->recv(eigenvals, 0, 1);
        }
        else
        {
            continue;
        }
    } 

    m_manager->widenotice(5) << "the size of my eigenvector matrix is: " << eigenvecs.rows() << "x" << eigenvecs.cols() << std::endl;
    m_manager->widenotice(5) << "the size of my eigenvalue vector is: " << eigenvals.rows() << "x" << eigenvals.cols() << std::endl;
    m_comm->barrier();
    //m_manager->widenotice(5) << "the size of my eigenvalue vector is: " << eigenvalue.rows() << "x" << eigenvalue.cols() << std::endl;
}

std::tuple<bool,double> SpectraHessian::getEigenvalue(unsigned int index)
{
    double eigenval = 0;
    div_t divres = div((int)nconv, (int)m_comm->getSizeGlobal());
    std::vector<int> counts(m_comm->getSizeGlobal(), divres.quot);
    for(int i = 0; i < m_comm->getSizeGlobal(); ++i)
    {
        if (i < divres.rem) 
            counts[i] += 1;
    }
    unsigned int starts = std::accumulate(counts.begin(), counts.begin()+m_comm->getRank(), 0);
    unsigned int ends = std::accumulate(counts.begin(), counts.begin()+m_comm->getRank()+1,0);
    if (index >= starts && index < ends)
    {
        eigenval = eigenvals[index-starts];
        return std::make_tuple(true,eigenval);
    }
    else
    {
        eigenval = 0;
        return std::make_tuple(false,eigenval);
    }
}

void SpectraHessian::saveEigenvalue(unsigned int index, std::shared_ptr<MPI::LogFile> logfile)
{
    std::tuple<bool, double> eigenval = getEigenvalue(index);
    if(std::get<0>(eigenval))
    {
        m_manager->notice(10) << "Process " << m_comm->getRank() << "is outputting eigenvalue" << std::endl;
        std::string outline = std::to_string(std::get<1>(eigenval)) + " ";
        logfile->write_shared(outline);
    }
}

std::tuple<bool,Eigen::VectorXd> SpectraHessian::getEigenvector(unsigned int index)
{
    Eigen::VectorXd eigenvector;
    div_t divres = div((int)nconv, (int)m_comm->getSizeGlobal());
    std::vector<int> counts(m_comm->getSizeGlobal(), divres.quot);
    for(int i = 0; i < m_comm->getSizeGlobal(); ++i)
    {
        if (i < divres.rem) 
            counts[i] += 1;
    }
    unsigned int starts = std::accumulate(counts.begin(), counts.begin()+m_comm->getRank(), 0);
    unsigned int ends = std::accumulate(counts.begin(), counts.begin()+m_comm->getRank()+1,0);//[sum(counts[:p+1]) for p in range(nprocs)]
    m_manager->widenotice(5) << starts << " " << index << " " << ends << std::endl;
    if (index >= starts && index < ends)
    {
        eigenvector = eigenvecs.col(index-starts);
        return std::make_tuple(true,eigenvector);
    }
    else
    {
        return std::make_tuple(false,eigenvector);
    }
} 

void SpectraHessian::saveEigenvector(unsigned int index, std::shared_ptr<MPI::LogFile> logfile)
{
    std::tuple<bool, Eigen::VectorXd> eigenvector = getEigenvector(index);
    if(std::get<0>(eigenvector))
    {
        m_manager->notice(10) << "Process " << m_comm->getRank() << "is outputting eigenvalue" << std::endl;
        int id_x = 0;
        for( auto p_i = m_sysdata->particles.begin(); p_i != m_sysdata->particles.end(); ++p_i)
        {
            //Input particle type 
            std::stringstream outline;
            outline << id_x << " "; 
            //Input particle position 
            outline << abr::get<position>(*p_i)[0] << " ";
            outline << abr::get<position>(*p_i)[1] << " ";
            outline << abr::get<position>(*p_i)[2] << " ";
            //Input the eigenvectors
            outline << std::get<1>(eigenvector)[2*id_x] << " ";
            outline << std::get<1>(eigenvector)[2*id_x+1] << " ";
            id_x += 1;
            //end the input
            outline << std::endl;
            logfile->write_shared(outline.str());
        }
    }
}
/*
void SpectraHessian::buildPseudoInverse()
{
    if (nconv > 0)
    {
        double cond = m_manager->pinv_tol;
        
        #pragma omp parallel for reduction(+:pseudoinverse) 
        for (unsigned int i = 0; i < nconv; ++i)
        {
            if (eigenvals[i] > cond)
            {
                pseudoinverse.noalias() += eigenvecs.col(i)*eigenvecs.col(i).transpose()/eigenvals[i];
            }
        }
    }
}
*/
void SpectraHessian::calculateNonAffineTensor() 
{
    if (eigenvals.size() > 0)
    {
        double cond = m_manager->pinv_tol;
        //#pragma omp parallel for reduction(+:nonaffinetensor) 
        for (unsigned int i = 0; i < eigenvals.size(); ++i)
        {
            if (eigenvals[i] > cond)
            {
                double laminv = 1/eigenvals[i];
                double xx = eigenvecs.col(i).dot(misforce.col(0));
                double yy = eigenvecs.col(i).dot(misforce.col(1));
                double xy = eigenvecs.col(i).dot(misforce.col(2));

                nonaffinetensor(0,0) += xx*xx*laminv; nonaffinetensor(0,1) += xx*yy*laminv; nonaffinetensor(0,2) += xx*xy*laminv;
                nonaffinetensor(1,0) += yy*xx*laminv; nonaffinetensor(1,1) += yy*yy*laminv; nonaffinetensor(1,2) += yy*xy*laminv;
                nonaffinetensor(2,0) += xy*xx*laminv; nonaffinetensor(2,1) += xy*yy*laminv; nonaffinetensor(2,2) += xy*xy*laminv;
            }
        }
        nonaffinetensor /= m_sysdata->simbox->vol;
        py::print(1,m_comm->getRank(),nonaffinetensor);
        std::vector< Eigen::MatrixXd > listofnonaffine(m_comm->getSizeGlobal());
        m_comm->gather_v(nonaffinetensor, listofnonaffine, 0);
        //MPI_Reduce(MPI_IN_PLACE,nonaffinetensor.data(),nonaffinetensor.rows(),MPI_DOUBLE,MPI_SUM,0,m_comm->getCommunicator());
        if (m_comm->isRoot())
        {
            nonaffinetensor = std::accumulate(listofnonaffine.begin(), listofnonaffine.end(), Eigen::MatrixXd::Zero(3,3).eval());
        }
        py::print(2,m_comm->getRank(),nonaffinetensor);
    }
}

std::tuple<bool,double> SpectraHessian::getNonAffineTensor(int i, int j) 
{
    if (m_comm->isRoot())
    {
        py::print(nonaffinetensor);
        return std::make_tuple(true,nonaffinetensor(i,j));
    }
    else
    {
        return std::make_tuple(false,0);
    }
}

void SpectraHessian::saveNonAffineTensor(unsigned int i, unsigned int j, std::shared_ptr<MPI::LogFile> logfile)
{
    std::tuple<bool, double> tensor_comp = getNonAffineTensor(i,j);
    if(std::get<0>(tensor_comp))
    {
        m_manager->notice(10) << "Process " << m_comm->getRank() << "is outputting nonaffinetensor" << std::endl;
        std::string outline = std::to_string(std::get<1>(tensor_comp)) + " ";
        logfile->write_shared(outline);
    }
}

void export_SpectraHessian(py::module& m)
{
    py::class_<SpectraHessian, std::shared_ptr<SpectraHessian> >(m,"SpectraHessian")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< HessianManager >, std::shared_ptr< MPI::Communicator > >())
    .def("setSystemData", &SpectraHessian::setSystemData)
    .def("buildHessianandMisForce", &SpectraHessian::buildHessianandMisForce)
    .def("getEigenPairs", &SpectraHessian::getEigenPairs)
    .def("getAllEigenPairs", &SpectraHessian::getAllEigenPairs)
    .def("getEigenvalue", &SpectraHessian::getEigenvalue)
    .def("saveEigenvalue", &SpectraHessian::saveEigenvalue)
    .def("getEigenvector", &SpectraHessian::getEigenvector)
    .def("saveEigenvector", &SpectraHessian::saveEigenvector)
    .def("calculateNonAffineTensor", &SpectraHessian::calculateNonAffineTensor)
    .def("getNonAffineTensor", &SpectraHessian::getNonAffineTensor)
    .def("saveNonAffineTensor", &SpectraHessian::saveNonAffineTensor)
    //.def("buildPseudoInverse", &SpectraHessian::buildPseudoInverse)
    ;
}

#endif //__SPECTRAHESSIAN_H__
