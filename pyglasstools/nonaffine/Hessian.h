#ifndef __HESSIAN_H__
#define __HESSIAN_H__

//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACKE

#include <vector>
#include <cmath>
#include <string>
#include <map>

#include <algorithm>
#include <memory>
#include <chrono>
#include <omp.h>

#include <pyglasstools/MathAndTypes.h>
#include <pyglasstools/ParticleSystem.h>
#include <pyglasstools/potential/PairPotential.h>

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
namespace spc = Spectra;

class PYBIND11_EXPORT Hessian
{
    public:
        SparseMatd hessian;
        Eigen::MatrixXd pseudoinverse;
        Eigen::VectorXd eigenvals;
        Eigen::MatrixXd eigenvecs;
        unsigned int nconv; 
        double frobeniuserror;    
        #pragma omp declare reduction( + : Eigen::MatrixXd : omp_out += omp_in ) initializer( omp_priv = Eigen::MatrixXd::Zero(omp_orig.rows(),omp_orig.cols()) )
        
        Hessian(std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential)
            : m_sysdata(sysdata), m_potential(potential)
            {
                double maxforce_rcut = potential->scaled_rcut;
                double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                                        std::end(abr::get<diameter>(m_sysdata->particles)) );
                maxforce_rcut *= maxdiameter;
                max_rcut = std::max(maxforce_rcut,maxdiameter); //whichever is maximum
                std::vector< Tripletd > hessian_triplet;
                hessian_length = (unsigned int)m_sysdata->simbox->dim*m_sysdata->particles.size();
                unsigned int estimate_nonzero_entries = (unsigned int)m_sysdata->simbox->dim*hessian_length*6;
                hessian_triplet.reserve(estimate_nonzero_entries);
                
                hessian.resize(hessian_length,hessian_length);
                pseudoinverse = Eigen::MatrixXd::Zero(hessian_length,hessian_length);
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

                            Eigen::Matrix3d offdiag_ij = -(m_potential->getBondStiffness()+m_potential->getPairForce())*nij*nij.transpose()
                                                         +Eigen::Matrix3d::Identity()*m_potential->getPairForce();
                            
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
                nconv = 0;
                frobeniuserror = 0;
            };
        virtual ~Hessian(){};

        void getFullEigenDecomposition(int maxiter, double tol)
        {
            // Construct matrix operation object using the wrapper class SparseGenMatProd
            spc::SparseSymShiftSolve<double> op(hessian);
            spc::SymEigsShiftSolver< double, spc::LARGEST_MAGN, spc::SparseSymShiftSolve<double> > eigs(&op, hessian_length-2, hessian_length, tol);
            eigs.init();
            py::print("Performing Full Eigendecomposition");
            nconv = eigs.compute(maxiter, tol,spc::SMALLEST_MAGN);
            int info = eigs.info();
            if(info == spc::SUCCESSFUL)
            {
                py::print("Eigenvalues and Eigenvectors computed successfully.");
                py::print("Number of converged eigenpairs is: ",nconv);
                py::print("Theoretical number of eigenpairs with non-zero eigenvalues is:",hessian_length-m_sysdata->simbox->dim);
                eigenvals = eigs.eigenvalues().real();
                eigenvecs = eigs.eigenvectors().real();
            }
        }   
        void buildPseudoInverse(double tol)
        {
            if (nconv > 0)
            {
                py::print("Assembling pseudoinverse . . .");
                if (nconv != hessian_length-m_sysdata->simbox->dim)
                {
                    py::print("[WARNING]! You are missing ",nconv-(hessian_length-m_sysdata->simbox->dim)," more converged eigenpairs");
                    py::print("[WARNING]! Check the norm errors at the end.");
                }
                double cond = hessian_length*eigenvals[nconv-1]*tol;
                #pragma omp parallel for reduction(+:pseudoinverse) 
                for (unsigned int i = 0; i < nconv; ++i)
                {
                    if (eigenvals[i] > cond)
                    {
                        pseudoinverse.noalias() += eigenvecs.col(i)*eigenvecs.col(i).transpose()/eigenvals[i];
                    }
                }
                py::print("Computing Error from ||AA^+A-A||");
                frobeniuserror = ((hessian*pseudoinverse-Eigen::MatrixXd::Identity(hessian_length,hessian_length))*hessian).norm();
            }
        } 
        void getEigenDecomposition(std::string selrule, int nev, int ncv, int maxiter, double tol)
        {
            // Construct matrix operation object using the wrapper class SparseGenMatProd
            int info = spc::NOT_COMPUTED;
            if (selrule == "LM")
            {
                spc::SparseGenMatProd<double> op(hessian);
                spc::SymEigsSolver< double, spc::LARGEST_MAGN, spc::SparseGenMatProd<double> > eigs(&op, nev, ncv);
                eigs.init();
                nconv = eigs.compute(maxiter, tol);
                info = eigs.info();
                eigenvals = eigs.eigenvalues().real();
                eigenvecs = eigs.eigenvectors().real();
            }    
            else if (selrule == "LR")
            {
                spc::SparseGenMatProd<double> op(hessian);
                spc::SymEigsSolver< double, spc::LARGEST_REAL, spc::SparseGenMatProd<double> > eigs(&op, nev, ncv);
                eigs.init();
                nconv = eigs.compute(maxiter, tol);
                info = eigs.info();
                eigenvals = eigs.eigenvalues().real();
                eigenvecs = eigs.eigenvectors().real();
            }    
            else if (selrule == "LA")
            {
                spc::SparseGenMatProd<double> op(hessian);
                spc::SymEigsSolver< double, spc::LARGEST_ALGE, spc::SparseGenMatProd<double> > eigs(&op, nev, ncv);
                eigs.init();
                nconv = eigs.compute(maxiter, tol);
                info = eigs.info();
                eigenvals = eigs.eigenvalues().real();
                eigenvecs = eigs.eigenvectors().real();
            }    
            else if (selrule == "SM")
            {
                spc::SparseSymShiftSolve<double> op(hessian);
                spc::SymEigsShiftSolver< double, spc::LARGEST_MAGN, spc::SparseSymShiftSolve<double> > eigs(&op, nev, ncv, tol);
                eigs.init();
                nconv = eigs.compute(maxiter, tol);
                info = eigs.info();
                eigenvals = eigs.eigenvalues().real();
                eigenvecs = eigs.eigenvectors().real();
            }    
            else if (selrule == "SR")
            {
                spc::SparseSymShiftSolve<double> op(hessian);
                spc::SymEigsShiftSolver< double, spc::LARGEST_REAL, spc::SparseSymShiftSolve<double> > eigs(&op, nev, ncv, tol);
                eigs.init();
                nconv = eigs.compute(maxiter, tol);
                info = eigs.info();
                eigenvals = eigs.eigenvalues().real();
                eigenvecs = eigs.eigenvectors().real();
            }    
            else if (selrule == "SA")
            {
                spc::SparseSymShiftSolve<double> op(hessian);
                spc::SymEigsShiftSolver< double, spc::LARGEST_ALGE, spc::SparseSymShiftSolve<double> > eigs(&op, nev, ncv, tol);
                eigs.init();
                nconv = eigs.compute();
                info = eigs.info();
                eigenvals = eigs.eigenvalues().real();
                eigenvecs = eigs.eigenvectors().real();
            }    
            else 
            {
                throw std::runtime_error("[ERROR] Selection Rule unrecognized");
            }

            if(info == spc::SUCCESSFUL)
            {
                py::print("Eigenvalues and Eigenvectors computed successfully.");
                py::print("Number of requested eigenpairs: ",nev);
                py::print("Number of converges eigenpairs: ",nconv);
            }
        }


    protected:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< PairPotential > m_potential;
        double max_rcut;
        unsigned int hessian_length;
};

void export_Hessian(py::module& m)
{
    py::class_<Hessian, std::shared_ptr<Hessian> >(m,"Hessian")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >  >())
    .def_readwrite("hessian", &Hessian::hessian)
    .def_readwrite("pseudoinverse", &Hessian::pseudoinverse)
    .def_readwrite("eigenvals", &Hessian::eigenvals)
    .def_readwrite("eigenvecs", &Hessian::eigenvecs)
    .def_readwrite("nconv", &Hessian::nconv)
    .def_readwrite("frobeniuserror", &Hessian::frobeniuserror)
    .def("getEigenDecomposition", &Hessian::getEigenDecomposition)
    .def("getFullEigenDecomposition", &Hessian::getFullEigenDecomposition)
    .def("buildPseudoInverse", &Hessian::buildPseudoInverse)
    ;
};

void export_SelectionRule(py::module& m)
{
    py::enum_<spc::SELECT_EIGENVALUE>(m, "SelectionRule")
        .value("LM", spc::SELECT_EIGENVALUE::LARGEST_MAGN) 	
        .value("LR", spc::SELECT_EIGENVALUE::LARGEST_REAL) 	
        .value("LI", spc::SELECT_EIGENVALUE::LARGEST_IMAG) 	
        .value("LA", spc::SELECT_EIGENVALUE::LARGEST_ALGE) 	
        .value("SM", spc::SELECT_EIGENVALUE::SMALLEST_MAGN) 	
        .value("SR", spc::SELECT_EIGENVALUE::SMALLEST_REAL) 	
        .value("SI", spc::SELECT_EIGENVALUE::SMALLEST_IMAG) 	
        .value("SA", spc::SELECT_EIGENVALUE::SMALLEST_ALGE) 	
        .value("BE", spc::SELECT_EIGENVALUE::BOTH_ENDS)
        .export_values();
}

/*
void export_SelectionRule(py::module& m)
{
    py::enum_<Pet::Kind>(pet, "Kind")
        .value("Dog", Pet::Kind::Dog)
        .value("Cat", Pet::Kind::Cat)
        .export_values();
}
*/
//Let's also export Spectra's enumeration types

#endif //__HESSIAN_H__
