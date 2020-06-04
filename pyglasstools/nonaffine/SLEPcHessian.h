#ifndef __SLEPCHESSIAN_H__
#define __SLEPCHESSIAN_H__

#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>
#include <memory>

#include <pyglasstools/MathAndTypes.h>
#include <pyglasstools/ParticleSystem.h>
#include <pyglasstools/potential/PairPotential.h>
#include "NonAffineManager.h"
#include "PETScObservable.h"

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
namespace py = pybind11;

#include <Aboria.h>
namespace abr = Aboria;

#include <slepceps.h>
namespace slepc
{
    PetscErrorCode initialize()
    {
        PetscErrorCode ierr = SlepcInitialize(NULL,NULL,(char*)0,NULL);
        return ierr;
    }

    void finalize()
    {
        SlepcFinalize();
    }
};

class PYBIND11_EXPORT SLEPcHessian
{
    protected:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< PairPotential > m_potential; //!< system pari potential
        std::shared_ptr< PETScManager > m_manager; //!< manager object, which handlers error/logging/etc.
        std::shared_ptr< MPI::Communicator > m_comm;  //!< MPI Communicator object, wrapping common MPI methods
        
        //List of observables to be computed!
        std::map<std::string, std::shared_ptr< PETScVectorFieldBase > > m_vectorfields;
        std::map<std::string, std::shared_ptr< PETScGlobalPropertyBase > > m_observables;
        
        double max_rcut;
        int nconv; 
        unsigned int hessian_length;
        double m_maxeigval;
       
        //PETSc objects for all computation 
        Mat hessian;
        EPS eps;
        Mat misforce;
        PetscErrorCode ierr;
        
        //Eigen::MatrixXd nonaffinetensor;
    public:
        double m_mineigval;

        SLEPcHessian(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                        std::shared_ptr< PETScManager > manager, std::shared_ptr< MPI::Communicator > comm);
        virtual ~SLEPcHessian()
        {
            MatDestroy(&hessian);
            MatDestroy(&misforce);
            EPSDestroy(&eps);
        };
        void setSystemData(std::shared_ptr<ParticleSystem> newsysdata)
        {
            m_sysdata = newsysdata;
        }
        
        void buildHessianandMisForce();
        
        Mat getHessian()
        {
            return hessian;
        } 
        
        virtual void addVectorField(const std::shared_ptr< PETScVectorFieldBase >& obs)
        {
            m_vectorfields.insert(std::pair<std::string, std::shared_ptr< PETScVectorFieldBase > >(obs->name, obs));
        }

        virtual void addGlobalProperty(const std::shared_ptr< PETScGlobalPropertyBase >& obs)
        {
            m_observables.insert(std::pair<std::string, std::shared_ptr< PETScGlobalPropertyBase > >(obs->name, obs));
        }
        //void getMaxEigenvalue_forMumps();
        void checkSmallestEigenvalue_forMumps();
        //void getAllEigenPairs_Mumps();
        //void getEigenPairs();
        //std::tuple<int,int> getPETScMatRange(unsigned int index); 
        
        //double getEigenvalue(unsigned int index); 
        //void saveEigenvalue(unsigned int index, std::shared_ptr<MPI::LogFile> logfile); 

        //std::vector<double> getEigenvector(unsigned int index, bool forall); 
        //void saveEigenvector(unsigned int index, std::shared_ptr<MPI::LogFile> logfile); 

        //void calculateNonAffineTensor(); 
        //void saveNonAffineTensor(unsigned int i, unsigned int j, std::shared_ptr<MPI::LogFile> logfile); 
       
        // Methods for force dipole calculations: 
        void solveForceDipoleProblem_Minimum();
        void solveForceDipoleProblem_Random();
        //void saveForceDipoleProblem(std::shared_ptr<MPI::LogFile> logfile);
};

SLEPcHessian::SLEPcHessian( std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                            std::shared_ptr<PETScManager> manager, std::shared_ptr<MPI::Communicator> comm)
    : m_sysdata(sysdata), m_potential(potential), m_manager(manager), m_comm(comm)
{
    int Dim = m_sysdata->simbox->dim;
    
    //Add command line options
    ierr = PetscOptionsInsertString(NULL,m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Constructing SLEPc Hessian Object \n");
    
    double maxforce_rcut = potential->scaled_rcut;
    double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                            std::end(abr::get<diameter>(m_sysdata->particles)) );
    maxforce_rcut *= maxdiameter;
    max_rcut = std::max(maxforce_rcut,maxdiameter); //whichever is maximum
    
    //Construct Hessian
    m_manager->printPetscNotice(5,"Assemble PETSc Sparse Matrix \n");
    nconv = 0;
    hessian_length = (unsigned int)Dim*m_sysdata->particles.size();
    
    buildHessianandMisForce();

    m_maxeigval = 0;
    m_mineigval = 0;
    //nonaffinetensor = Eigen::MatrixXd::Zero((unsigned int)Dim*((unsigned int)Dim+1)/2,(unsigned int)Dim*((unsigned int)Dim+1)/2);                
};

void SLEPcHessian::buildHessianandMisForce()
{
        int Dim = m_sysdata->simbox->dim;

        //Construct the hessian matrix!
        MatCreate(PETSC_COMM_WORLD,&hessian);
        MatSetSizes(hessian,PETSC_DETERMINE,PETSC_DETERMINE,hessian_length,hessian_length);
        MatSetType(hessian,MATAIJ);
        MatSetUp(hessian);

        //Construct mismatch force vector
        MatCreate(PETSC_COMM_WORLD,&misforce);
        MatSetSizes(misforce,PETSC_DETERMINE,PETSC_DETERMINE,hessian_length,(unsigned int)Dim*((unsigned int)Dim+1)/2);
        MatSetType(misforce,MATDENSE);
        MatSetUp(misforce);

        PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(hessian, &Istart, &Iend);
        for( auto p_i = m_sysdata->particles.begin(); p_i != m_sysdata->particles.end(); ++p_i)
        {
            //py::print(id_i,id_i+1,Iend);
            for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                        abr::get<position>(*p_i), max_rcut); p_j != false; ++p_j)
            {
                PetscInt id_j = abr::get<abr::id>(*p_j);
                PetscInt id_i = abr::get<abr::id>(*p_i);
                
                if (id_i != id_j && 2*id_i < Iend && 2*id_i >= Istart && 2*id_i + 1 < Iend && 2*id_i+1 >= Istart)
                {
                    //Compute a list of local obsercavles 
                    //Next we loop through the j-th particles for the virial stress
                    //set the distance between particle i and particle j
                    Eigen::Vector3d rij(-p_j.dx()[0], -p_j.dx()[1], -p_j.dx()[2]);
                    Eigen::Vector3d nij = rij.normalized();
                    
                    double di =  abr::get<diameter>(*p_i);
                    double dj =  abr::get<diameter>(*p_j);
                    //Don't forget to set diameters of the potential
                    
                       
                    //Make sure the particle is unique
                    if (m_potential->getRcut(rij, di, dj) > rij.dot(rij)) 
                    {
                        double factor = m_potential->getBondStiffness(rij, di, dj)+m_potential->getPairForce(rij, di, dj);
                        Eigen::Matrix3d offdiag_ij = -factor*nij*nij.transpose()
                                                     +Eigen::Matrix3d::Identity()*m_potential->getPairForce(rij, di, dj);

                        MatSetValue(misforce,2*id_i,0, -factor*rij[0]*nij[0]*nij[0],ADD_VALUES); 
                        MatSetValue(misforce,2*id_i+1,0, -factor*rij[0]*nij[0]*nij[1],ADD_VALUES); 
                        
                        MatSetValue(misforce,2*id_i,1, -factor*rij[1]*nij[1]*nij[0],ADD_VALUES); 
                        MatSetValue(misforce,2*id_i+1,1, -factor*rij[1]*nij[1]*nij[1],ADD_VALUES); 
                        
                        MatSetValue(misforce,2*id_i,2, -factor*rij[0]*nij[1]*nij[0],ADD_VALUES); 
                        MatSetValue(misforce,2*id_i+1,2, -factor*rij[0]*nij[1]*nij[1],ADD_VALUES); 
                        
                        MatSetValue(hessian,2*id_i, 2*id_j, offdiag_ij(0,0),ADD_VALUES);
                        MatSetValue(hessian,2*id_i, 2*id_j+1, offdiag_ij(0,1) ,ADD_VALUES);
                        MatSetValue(hessian,2*id_i+1, 2*id_j, offdiag_ij(1,0) ,ADD_VALUES);
                        MatSetValue(hessian,2*id_i+1, 2*id_j+1, offdiag_ij(1,1) ,ADD_VALUES);
                        
                        MatSetValue(hessian,2*id_i, 2*id_i, -offdiag_ij(0,0) ,ADD_VALUES);
                        MatSetValue(hessian,2*id_i, 2*id_i+1, -offdiag_ij(0,1) ,ADD_VALUES);
                        MatSetValue(hessian,2*id_i+1, 2*id_i, -offdiag_ij(1,0) ,ADD_VALUES);
                        MatSetValue(hessian,2*id_i+1, 2*id_i+1, -offdiag_ij(1,1) ,ADD_VALUES);
                    }
                }
            }
        }

        MatAssemblyBegin(hessian,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(hessian,MAT_FINAL_ASSEMBLY);
        
        MatAssemblyBegin(misforce,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(misforce,MAT_FINAL_ASSEMBLY);
        
        EPSCreate(PETSC_COMM_WORLD,&eps);
        ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        EPSSetOperators(eps,hessian,NULL);
};

void SLEPcHessian::solveForceDipoleProblem_Minimum()
{
    if (m_mineigval >= 0)
    {
        Vec forcing; 
        KSP ksp;
        m_vectorfields["forcedipole"]->createVector(hessian);
        MatCreateVecs(hessian,NULL,&forcing);
        PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(hessian, &Istart, &Iend);
        
        //choose a random particle
        double forcedipole = std::numeric_limits<double>::max();
        PetscInt id_j, id_i;
        Eigen::Vector3d rij; 
        //All processes will look for
        m_manager->notice(5) << "Looking for pairs of particles to perturb" << std::endl; 
        for( auto p_i = m_sysdata->particles.begin(); p_i != m_sysdata->particles.end(); ++p_i)
        {
            //std::cout << dist6(rng) << std::endl;
            //py::print(id_i,id_i+1,Iend);
            for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                        abr::get<position>(*p_i), max_rcut); p_j != false; ++p_j)
            {
                PetscInt tempid_j = abr::get<abr::id>(*p_j);
                PetscInt tempid_i = abr::get<abr::id>(*p_i);
               
                //This needs to be changed so that it works in 2/3 dimensions 
                if (tempid_i < tempid_j)
                {
                    //Compute a list of local obsercavles 
                    //Next we loop through the j-th particles for the virial stress
                    //set the distance between particle i and particle j
                    Eigen::Vector3d bar_rij(-p_j.dx()[0], -p_j.dx()[1], -p_j.dx()[2]);
                    
                    //Don't forget to set diameters of the potential
                    double di =  abr::get<diameter>(*p_i);
                    double dj =  abr::get<diameter>(*p_j);
                    
                    //Make sure the particle is unique
                    //if (m_potential->getRcut() > rij.dot(rij) && abs(m_potential->getPairForce()) < abs(force_threshold))
                    double testdipole = abs(m_potential->getPairForce(bar_rij,di,dj));
                    if (m_potential->getRcut(bar_rij,di,dj) > bar_rij.dot(bar_rij) && abs(forcedipole) > testdipole && testdipole > 0 )
                    {
                        forcedipole = m_potential->getPairForce(bar_rij,di,dj);
                        rij = bar_rij;
                        id_i = tempid_i; 
                        id_j = tempid_j; 
                    }
                }
            }
        }

        //Now determine which rank has the smallest force. Once determined, we will broadcast the values
        m_manager->notice(5) << "Assembling forcing vector" << std::endl;
        std::stringstream streamm;
        streamm << "Perturbing particle " << id_i << " and particle " << id_j << " with force_x " << forcedipole*rij[0] << " with force_y " << forcedipole*rij[1] << std::endl; 
        m_manager->printPetscNotice(5, streamm.str());
        if (id_i != id_j && 2*id_i < Iend && 2*id_i >= Istart && 2*id_i + 1 < Iend && 2*id_i+1 >= Istart)
        {
            VecSetValue(forcing,2*id_i,-forcedipole*rij[0], ADD_VALUES);
            VecSetValue(forcing,2*id_i+1,-forcedipole*rij[1], ADD_VALUES);
            VecSetValue(forcing,2*id_j,forcedipole*rij[0], ADD_VALUES);
            VecSetValue(forcing,2*id_j+1,forcedipole*rij[1], ADD_VALUES);
        }
        VecAssemblyBegin(forcing);
        VecAssemblyEnd(forcing);
        
        Eigen::Vector3d center;
        center <<  0.5*(2*abr::get<position>(m_sysdata->particles[id_i])[0]-rij[0]),
                   0.5*(2*abr::get<position>(m_sysdata->particles[id_i])[1]-rij[1]),
                   0.5*(2*abr::get<position>(m_sysdata->particles[id_i])[2]-rij[2]);
        for (int i = 0; i < center.size(); ++i)
        {
            double x_size = m_sysdata->simbox->boxsize[i];
            if (center(i) <  -x_size * 0.5)
                center(i) += x_size;
            else if (center(i) >=  x_size * 0.5) 
                center(i) -= x_size;
            m_observables["forcedipole_center"]->setValue(center[i],i);
            m_observables["forcedipole_pi"]->setValue(abr::get<position>(m_sysdata->particles[id_i])[i],i);
            m_observables["forcedipole_pj"]->setValue(abr::get<position>(m_sysdata->particles[id_j])[i],i);
        }
        //set the nullspace
        m_manager->printPetscNotice(5, "Setting the Null Space of our Hessian. \n");
        MatNullSpace constant;
        MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&constant);
        MatSetNullSpace(hessian,constant);

        //set the values according to particleid
        KSPCreate(PETSC_COMM_WORLD,&ksp);
        KSPSetOperators(ksp,hessian,hessian);
        KSPSetFromOptions(ksp);
        m_manager->printPetscNotice(5, "Solve the force dipole problem. \n");
        KSPSolve(ksp,forcing,m_vectorfields["forcedipole"]->vectorobs);
        MatNullSpaceDestroy(&constant);
        m_manager->printPetscNotice(5, "Force dipole calculation finished. \n");
    }
};

void SLEPcHessian::solveForceDipoleProblem_Random()
{
    if (m_mineigval >= 0)
    {
        Vec forcing; 
        KSP ksp;
        m_vectorfields["forcedipole"]->createVector(hessian);
        MatCreateVecs(hessian,NULL,&forcing);
        PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(hessian, &Istart, &Iend);
        
        //choose a random particle
        double forcedipole = std::numeric_limits<double>::max();
        PetscInt id_j, id_i;
        Eigen::Vector3d rij; 
        //All processes will look for
        if (m_comm->isRoot())
        {
            std::random_device dev;
            std::mt19937 rng(dev());
            std::uniform_int_distribution<std::mt19937::result_type> dist6(0,m_sysdata->particles.size()-1);
            m_manager->notice(5) << "Looking for pairs of particles to perturb" << std::endl; 
            //for( auto p_i = m_sysdata->particles.begin(); p_i != m_sysdata->particles.end(); ++p_i)
            bool notfound = true;
            do
            {
                //std::cout << dist6(rng) << std::endl;
                //py::print(id_i,id_i+1,Iend);
                auto p_i = m_sysdata->particles.begin()+dist6(rng);
                for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                            abr::get<position>(*p_i), max_rcut); p_j != false; ++p_j)
                {
                    PetscInt tempid_j = abr::get<abr::id>(*p_j);
                    PetscInt tempid_i = abr::get<abr::id>(*p_i);
                   
                    //This needs to be changed so that it works in 2/3 dimensions 
                    if (tempid_i < tempid_j)
                    {
                        //Compute a list of local obsercavles 
                        //Next we loop through the j-th particles for the virial stress
                        //set the distance between particle i and particle j
                        Eigen::Vector3d bar_rij(-p_j.dx()[0], -p_j.dx()[1], -p_j.dx()[2]); //this is rij = ri-rj, with PBC taken into account
                        
                        //Don't forget to set diameters of the potential
                        double di =  abr::get<diameter>(*p_i);
                        double dj =  abr::get<diameter>(*p_j);
                        
                        //Make sure the particle is unique
                        //if (m_potential->getRcut() > rij.dot(rij) && abs(m_potential->getPairForce()) < abs(force_threshold))
                        double testdipole = abs(m_potential->getPairForce(bar_rij,di,dj));
                        if (    m_potential->getRcut(bar_rij,di,dj) > bar_rij.dot(bar_rij) && 
                                abs(forcedipole) > testdipole && testdipole*bar_rij.norm() > m_manager->fd_random_min 
                                && testdipole*bar_rij.norm() < m_manager->fd_random_max)
                        {
                            m_manager->notice(5) << "Found it!" << std::endl;
                            m_manager->notice(5) << "Force min: " << m_manager->fd_random_min << " and Force max: " << m_manager->fd_random_max << std::endl;
                            forcedipole = m_potential->getPairForce(bar_rij,di,dj);
                            rij = bar_rij;
                            id_i = tempid_i; 
                            id_j = tempid_j;
                            notfound = false; 
                            break;
                        }
                    }
                }
            }
            while (notfound);
        }
        m_comm->barrier();
        //Now, broadcast the results for id_i, id_j, rij, and forcedipole
        //Now determine which rank has the smallest force. Once determined, we will broadcast the values
        m_comm->bcast(id_i,0); 
        m_comm->bcast(id_j,0); 
        m_comm->bcast(forcedipole,0); 
        m_comm->bcast(rij,0);

        m_manager->notice(5) << "Assembling forcing vector" << std::endl;
        std::stringstream streamm;
        streamm << "Perturbing particle " << id_i << " and particle " << id_j << " with force_x " << forcedipole*rij[0] << " with force_y " << forcedipole*rij[1] << std::endl; 
        m_manager->printPetscNotice(5, streamm.str());
        if (id_i != id_j && 2*id_i < Iend && 2*id_i >= Istart && 2*id_i + 1 < Iend && 2*id_i+1 >= Istart)
        {
            VecSetValue(forcing,2*id_i,-forcedipole*rij[0], ADD_VALUES);
            VecSetValue(forcing,2*id_i+1,-forcedipole*rij[1], ADD_VALUES);
            VecSetValue(forcing,2*id_j,forcedipole*rij[0], ADD_VALUES);
            VecSetValue(forcing,2*id_j+1,forcedipole*rij[1], ADD_VALUES);
        }
        VecAssemblyBegin(forcing);
        VecAssemblyEnd(forcing);
        
        Eigen::Vector3d center;
        center <<  0.5*(2*abr::get<position>(m_sysdata->particles[id_i])[0]-rij[0]),
                   0.5*(2*abr::get<position>(m_sysdata->particles[id_i])[1]-rij[1]),
                   0.5*(2*abr::get<position>(m_sysdata->particles[id_i])[2]-rij[2]);
        for (int i = 0; i < center.size(); ++i)
        {
            double x_size = m_sysdata->simbox->boxsize[i];
            if (center(i) <  -x_size * 0.5)
                center(i) += x_size;
            else if (center(i) >=  x_size * 0.5) 
                center(i) -= x_size;
            m_observables["forcedipole_center"]->setValue(center[i],i);
            m_observables["forcedipole_pi"]->setValue(abr::get<position>(m_sysdata->particles[id_i])[i],i);
            m_observables["forcedipole_pj"]->setValue(abr::get<position>(m_sysdata->particles[id_j])[i],i);
        }

        //set the nullspace
        m_manager->printPetscNotice(5, "Setting the Null Space of our Hessian. \n");
        MatNullSpace constant;
        Vec nullvec[2];//, nully
        MatCreateVecs(hessian,NULL,&nullvec[0]);
        MatCreateVecs(hessian,NULL,&nullvec[1]);
        for (int i = 0; i < m_sysdata->particles.size()*m_sysdata->simbox->dim; ++i)
        {
            VecSetValue(nullvec[0],2*id_i,1, INSERT_VALUES);
            VecSetValue(nullvec[0],2*id_i+1,0, INSERT_VALUES);
            VecSetValue(nullvec[1],2*id_i,0, INSERT_VALUES);
            VecSetValue(nullvec[1],2*id_i+1,1, INSERT_VALUES);
        }
        for (int i = 0; i < m_sysdata->simbox->dim; ++i)
        {
            VecAssemblyBegin(nullvec[i]);
            VecAssemblyEnd(nullvec[i]);
            VecNormalize(nullvec[i],NULL);
        }
        MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_FALSE,2,nullvec,&constant);
        MatSetNullSpace(hessian,constant);

        //set the values according to particleid
        KSPCreate(PETSC_COMM_WORLD,&ksp);
        KSPSetOperators(ksp,hessian,hessian);
        KSPSetFromOptions(ksp);
        m_manager->printPetscNotice(5, "Solve the force dipole problem. \n");
        KSPSolve(ksp,forcing,m_vectorfields["forcedipole"]->vectorobs);
        MatNullSpaceDestroy(&constant);
        m_manager->printPetscNotice(5, "Force dipole calculation finished. \n");
    }
};

void SLEPcHessian::checkSmallestEigenvalue_forMumps()
{
    //This should only be done when computing a small amount of eigenvalues
    m_manager->printPetscNotice(5,"Computing Smallest Eigenvalue \n");
    //and setting the which to anything other than EPS_ALL
    EPSSetDimensions(eps,1,m_sysdata->particles.size(),m_sysdata->particles.size());
    EPSSetTolerances(eps,PETSC_DEFAULT,PETSC_DEFAULT); 
    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    EPSSetUp(eps);
    EPSSolve(eps);
    
    PetscReal vr;
    ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSGetEigenvalue(eps,0,&vr,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Smallest Eigenvalue: "+detail::to_string_sci(vr)+"\n");
    m_mineigval = vr;
    if (m_manager->getNoticeLevel() >= 6)
    {
        m_manager->printPetscNotice(6,"Summary of Normal Mode Analysis: \n");
        ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
};
/*
template<int Dim>
double SLEPcHessian<Dim>::getEigenvalue(unsigned int index) 
{
        PetscReal lambda_r;//Vec xr,global_xr;
        ierr = EPSGetEigenvalue(eps,index,&lambda_r,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        return lambda_r;
};

template<int Dim>
void SLEPcHessian<Dim>::saveEigenvalue(unsigned int index, std::shared_ptr<MPI::LogFile> logfile)
{
    double eigenval = getEigenvalue(index);
    if(m_comm->isRoot())
    {
        m_manager->notice(10) << "Process " << m_comm->getRank() << "is outputting eigenvalue" << std::endl;
        std::string outline = std::to_string(eigenval) + " ";
        logfile->write_shared(outline);
    }
}

template<int Dim>
std::vector<double> SLEPcHessian<Dim>::getEigenvector(unsigned int index, bool forall) 
{
        int globsize;
        Vec xr,global_xr;
        VecScatter ctx;
        MatCreateVecs(hessian,NULL,&xr);
        VecGetSize(xr,&globsize);
        std::vector<double> eigenvector;
        if (forall)
            VecScatterCreateToAll(xr,&ctx,&global_xr);
        else 
            VecScatterCreateToZero(xr,&ctx,&global_xr);
        
        double *tempvecp = new double[globsize];
        VecGetArray(global_xr,&tempvecp);
        ierr = EPSGetEigenvector(eps,index,xr,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        VecScatterBegin(ctx,xr,global_xr,INSERT_VALUES,SCATTER_FORWARD);
        VecScatterEnd(ctx,xr,global_xr,INSERT_VALUES,SCATTER_FORWARD);
        if (forall)
        {
            std::copy(tempvecp, tempvecp+globsize, std::back_inserter(eigenvector));
        } 
        else
        {
            if (m_comm->isRoot())
            {
                std::copy(tempvecp, tempvecp+globsize, std::back_inserter(eigenvector));
            }
            else
            {
                eigenvector.resize(globsize);
            }
        }
        //Destroy all objects
        VecRestoreArray(global_xr,&tempvecp);
        delete [] tempvecp;
        VecScatterDestroy(&ctx);
        VecDestroy(&global_xr);
        VecDestroy(&xr);
        return eigenvector;
};

template<int Dim>
void SLEPcHessian<Dim>::saveEigenvector(unsigned int index, std::shared_ptr<MPI::LogFile> logfile)
{
    std::vector<double> eigenvector = getEigenvector(index,false);
    m_manager->notice(10) << "Process " << m_comm->getRank() << "is outputting eigenvector" << std::endl;
    if (m_comm->isRoot())
    {
        int id_x = 0;
        m_manager->notice(10) << "didwereachhere? " << m_comm->getRank() << "is outputting eigenvalue" << std::endl;
        for( auto p_i = m_sysdata->particles.begin(); p_i != m_sysdata->particles.end(); ++p_i)
        {
            //Input particle type 
            std::stringstream outline;
            outline << id_x << " "; 
            //Input particle position 
            outline << abr::get<position>(m_sysdata->particles[id_x])[0] << " ";
            outline << abr::get<position>(m_sysdata->particles[id_x])[1] << " ";
            outline << abr::get<position>(m_sysdata->particles[id_x])[2] << " ";
            //Input the eigenvectors
            outline << eigenvector[2*id_x] << " ";
            outline << eigenvector[2*id_x+1] << " ";
            //end the input
            outline << std::endl;
            id_x+=1;
            logfile->write_shared(outline.str());
        }
    }
    m_comm->barrier();
}

template<int Dim>
std::tuple<int,int> SLEPcHessian<Dim>::getPETScMatRange(unsigned int index) 
{
        Vec xr;
        MatCreateVecs(hessian,NULL,&xr);
        int Istart,Iend;
        VecGetOwnershipRange(xr,&Istart, &Iend);
        VecDestroy(&xr);
        return std::make_tuple(Istart,Iend);
};

template<int Dim>
void SLEPcHessian<Dim>::calculateNonAffineTensor() 
{
    if (nconv > 0)
    {
        ierr = EPSGetEigenvalue(eps,nconv-1,&m_maxeigval,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        double cond = m_manager->pinv_tol;//tol;
        Vec xr,mult,seqmult;
        VecScatter ctx;
        
        PetscReal kr,re;
        MatCreateVecs(hessian,NULL,&xr);
        MatCreateVecs(misforce,&mult,NULL);
        //Create A Scatterrer and the resulting sequential vector holding temporary values 
        VecScatterCreateToAll(mult,&ctx,&seqmult);
        
        double *tempvecp = new double[nonaffinetensor.rows()];//tempvec.data();
        VecGetArray(seqmult,&tempvecp);

        for (int i = 0; i < nconv; ++i)
        {   
            ierr = EPSGetEigenpair(eps,i,&kr,NULL,xr,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            #if defined(PETSC_USE_COMPLEX)
                re = PetscRealPart(kr);
            #else
                re = kr;
            #endif
            //Get the info here
            if (re > cond)
            {
                double laminv = 1/re;
                ierr = MatMultTranspose(misforce,xr,mult); CHKERRABORT(PETSC_COMM_WORLD,ierr);
                
                VecScatterBegin(ctx,mult,seqmult,INSERT_VALUES,SCATTER_FORWARD);
                VecScatterEnd(ctx,mult,seqmult,INSERT_VALUES,SCATTER_FORWARD);
                
                nonaffinetensor(0,0) += tempvecp[0]*tempvecp[0]*laminv; 
                nonaffinetensor(0,1) += tempvecp[0]*tempvecp[1]*laminv; 
                nonaffinetensor(0,2) += tempvecp[0]*tempvecp[2]*laminv;

                nonaffinetensor(1,0) += tempvecp[1]*tempvecp[0]*laminv; 
                nonaffinetensor(1,1) += tempvecp[1]*tempvecp[1]*laminv; 
                nonaffinetensor(1,2) += tempvecp[1]*tempvecp[2]*laminv;
                
                nonaffinetensor(2,0) += tempvecp[2]*tempvecp[0]*laminv; 
                nonaffinetensor(2,1) += tempvecp[2]*tempvecp[1]*laminv; 
                nonaffinetensor(2,2) += tempvecp[2]*tempvecp[2]*laminv;
            }
        }
        VecRestoreArray(seqmult,&tempvecp);
        delete [] tempvecp;
        VecScatterDestroy(&ctx);
        VecDestroy(&seqmult);
        VecDestroy(&mult);
        nonaffinetensor /= m_sysdata->simbox->vol;
    }
};

template<int Dim>
void SLEPcHessian<Dim>::saveNonAffineTensor(unsigned int i, unsigned int j, std::shared_ptr<MPI::LogFile> logfile)
{
    if(m_comm->isRoot())
    {
        m_manager->notice(10) << "Process " << m_comm->getRank() << "is outputting nonaffinetensor" << std::endl;
        std::string outline = std::to_string(nonaffinetensor(i,j)) + " ";
        logfile->write_shared(outline);
    }
}

template<int Dim>
void SLEPcHessian<Dim>::getAllEigenPairs_Mumps()
{
    getMaxEigenvalue_forMumps();

    //Set value MUMPS work-space
    ierr = PetscOptionsInsertString(NULL,"-mat_mumps_icntl_14 200");CHKERRABORT(PETSC_COMM_WORLD,ierr);
    // Set solver parameters at runtime
    ierr = EPSSetFromOptions(eps);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    
    EPSType type;
    ST             st;           
    KSP            ksp;
    PC             pc;
    PetscInt       maxit;
    PetscReal      inta,intb,tol;
    
    ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    
    //    Spectrum slicing requires Krylov-Schur
    ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSGetType(eps,&type);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Solution method: "+std::string(type)+"\n");
    
    ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Stopping condition: tol="+detail::to_string_sci(tol)+", maxit="+std::to_string(maxit)+"\n");
    //
    //   Set interval for spectrum slicing
    //
    inta = m_manager->lowerbound_tol; //tol;
    intb = m_maxeigval*(1+inta);//PETSC_MAX_REAL;
    ierr = EPSSetInterval(eps,inta,intb);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Search interval is: ["+detail::to_string_sci(inta)+", "+ detail::to_string_sci(intb)+"]\n");
    
    //
    //     Set shift-and-invert with Cholesky; select MUMPS if available
    //
    ierr = EPSGetST(eps,&st);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = STSetType(st,STSINVERT);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = STGetKSP(st,&ksp);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = KSPSetType(ksp,KSPPREONLY);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = PCSetType(pc,PCCHOLESKY);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    #if defined(PETSC_HAVE_MUMPS)
        #if defined(PETSC_USE_COMPLEX)
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Spectrum slicing with MUMPS is not available for complex scalars");
        #endif
        ierr = EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE);CHKERRABORT(PETSC_COMM_WORLD,ierr);  // enforce zero detection 
        ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        //
        //Add several MUMPS options (see ex43.c for a better way of setting them in program):
        //'-mat_mumps_icntl_13 1': turn off ScaLAPACK for matrix inertia
        //'-mat_mumps_icntl_24 1': detect null pivots in factorization (for the case that a shift is equal to an eigenvalue)
        //'-mat_mumps_cntl_3 <tol>': a tolerance used for null pivot detection (must be larger than machine epsilon)

        //Note: depending on the interval, it may be necessary also to increase the workspace:
        //'-mat_mumps_icntl_14 <percentage>': increase workspace with a percentage (50, 100 or more)
        //
        ierr = PetscOptionsInsertString(NULL,"-mat_mumps_icntl_13 1");CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscOptionsInsertString(NULL,"-mat_mumps_icntl_24 1");CHKERRABORT(PETSC_COMM_WORLD,ierr);
        std::string cntl_3 = "-mat_mumps_cntl_3 ";
        cntl_3 +=  std::to_string(std::numeric_limits<float>::epsilon()*inta);
        ierr = PetscOptionsInsertString(NULL,cntl_3.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
    #endif

    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Solve the eigensystem
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ierr = EPSSetUp(eps);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Solving for Eigenpairs . . . \n");
    ierr = EPSSolve(eps); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    
    ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                Display solution and clean up
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    //Get number of converged approximate eigenpairs
    ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Number of converged eigenpairs: "+std::to_string(nconv)+"\n");
    if (m_manager->getNoticeLevel() >= 6)
    {
        m_manager->printPetscNotice(6,"Summary of Normal Mode Analysis: \n");
        ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
};

template<int Dim>
void SLEPcHessian<Dim>::getMaxEigenvalue_forMumps()
{
    //This should only be done when computing a small amount of eigenvalues
    m_manager->printPetscNotice(5,"Computing Maximum Eigenvalue \n");
    //and setting the which to anything other than EPS_ALL
    EPSSetDimensions(eps,1,m_sysdata->particles.size(),m_sysdata->particles.size());
    EPSSetTolerances(eps,PETSC_DEFAULT,PETSC_DEFAULT); 
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
    EPSSetUp(eps);
    EPSSolve(eps);
    
    ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSGetEigenvalue(eps,0,&m_maxeigval,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Maximum Eigenvalue: "+detail::to_string_sci(m_maxeigval)+"\n");
};


template<int Dim>
void SLEPcHessian<Dim>::getEigenPairs()
{
    //This should only be done when computing a small amount of eigenvalues
    //and setting the which to anything other than EPS_ALL
    EPSType type;
    Vec            xr,xi;
    PetscReal      tol;
    PetscInt       nev,maxit,its;
    EPSWhich whicheig;
    
    EPSSetFromOptions(eps);
    EPSGetWhichEigenpairs(eps, &whicheig);
    if (whicheig == EPS_ALL)
    {
        m_manager->printPetscNotice(5,"Performing spectrum slicing! Make sure all options are configured correctly \n");
    }
    
    EPSSetUp(eps);
    EPSSolve(eps);
    //Optional: Get some information from the solver and display it
    ierr = EPSGetIterationNumber(eps,&its);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Number of iterations of the method: "+std::to_string(its)+"\n");
    ierr = EPSGetType(eps,&type);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Solution method: "+std::string(type)+"\n");
    ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Number of requested eigenvalues: "+std::to_string(nev)+"\n");
    ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Stopping condition: tol="+detail::to_string_sci(tol)+", maxit="+std::to_string(maxit)+"\n");
    
    ierr = MatCreateVecs(hessian,NULL,&xr);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = MatCreateVecs(hessian,NULL,&xi);CHKERRABORT(PETSC_COMM_WORLD,ierr);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                Display solution and clean up
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    //
    // Get number of converged approximate eigenpairs
    //
    ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Number of converged eigenpairs: "+std::to_string(nconv)+"\n");
    if (m_manager->getNoticeLevel() >= 6)
    {
        m_manager->printPetscNotice(6,"Summary of Normal Mode Analysis: \n");
        ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
    ierr = VecDestroy(&xr);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecDestroy(&xi);CHKERRABORT(PETSC_COMM_WORLD,ierr);
};
*/

void export_SLEPcHessian(py::module& m)
{
    py::class_<SLEPcHessian, std::shared_ptr<SLEPcHessian> >(m,"SLEPcHessian")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< PETScManager>  , std::shared_ptr< MPI::Communicator > >())
    .def("setSystemData", &SLEPcHessian::setSystemData)
    .def("buildHessianandMisForce", &SLEPcHessian::buildHessianandMisForce)
    .def("solveForceDipoleProblem_Minimum", &SLEPcHessian::solveForceDipoleProblem_Minimum)
    .def("solveForceDipoleProblem_Random", &SLEPcHessian::solveForceDipoleProblem_Random)
    .def("checkSmallestEigenvalue_forMumps",&SLEPcHessian::checkSmallestEigenvalue_forMumps)
    .def("addVectorField", &SLEPcHessian::addVectorField)
    .def("addGlobalProperty", &SLEPcHessian::addGlobalProperty)
    //.def("saveForceDipoleProblem", &T::saveForceDipoleProblem)
    //.def("getEigenPairs", &T::getEigenPairs)
    //.def("getAllEigenPairs", &T::getAllEigenPairs_Mumps)
    //.def("calculateNonAffineTensor", &T::calculateNonAffineTensor)
    //.def("saveNonAffineTensor", &T::saveNonAffineTensor)
    //.def("saveEigenvector", &T::saveEigenvector)
    //.def("getEigenvalue", &T::getEigenvalue)
    //.def("saveEigenvalue", &T::saveEigenvalue)
    //.def("getPETScMatRange", &T::getPETScMatRange)
    .def_readwrite("mineigval", &SLEPcHessian::m_mineigval)
    ;
};



//CHKERRABORT(PETSC_COMM_WORLD,ierr);
//Getits global size and divide by dimensionality
//VecGetSize(vectorobs,&globsize);
//Now, let's compute the allocation for the split vector
/*
div_t divres = div((int)globsize, (int)m_comm->getSizeGlobal());
std::vector<int> counts(m_comm->getSizeGlobal(), divres.quot);

for(unsigned int i = 0; i < m_comm->getSizeGlobal(); ++i)
{
    if (i < (unsigned int)divres.rem) 
        counts[i] += 1;
}
PetscInt starts = std::accumulate(counts.begin(), counts.begin()+m_comm->getRank(), 0);
ISCreateStride(m_comm->getCommunicator(),counts[m_comm->getRank()],starts,Dim,&is[0]);
ISCreateStride(m_comm->getCommunicator(),counts[m_comm->getRank()],starts+1,Dim,&is[1]);
if (Dim == 3)
    ISCreateStride(m_comm->getCommunicator(),counts[m_comm->getRank()],starts+Dim-1,Dim,&is[Dim-1]);
*/
/*
for(unsigned int i = 0; i < Dim; ++i)
{
    VecRestoreSubVector(vectorobs,is[i],&split_vectorobs[i]);
    VecDestroy(&split_vectorobs[i]);
    ISDestroy(&is[i]);
}
*/
/*
virtual void splitVector()
{
VecGetSubVector(vectorobs,is[0],&split_vectorobs[0]);
VecGetSubVector(vectorobs,is[1],&split_vectorobs[1]);
if (Dim == 3)
   VecGetSubVector(vectorobs,is[Dim-1],&split_vectorobs[Dim-1]);
}
*/
#endif //__SLEPCHESSIAN_H__
