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
        std::shared_ptr< PairPotential > m_potential;
        std::shared_ptr< PETScManager > m_manager;
        std::shared_ptr< MPI::Communicator > m_comm; 
        double max_rcut;
        unsigned int hessian_length;
        double m_maxeigval;
        
        Mat hessian;
        EPS eps;
        Mat misforce;
        PetscErrorCode ierr;
        
        Eigen::MatrixXd nonaffinetensor;
        Vec forcedipole_disp;
    public:
        int nconv; 

        SLEPcHessian(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                        std::shared_ptr<PETScManager> manager, std::shared_ptr< MPI::Communicator > comm);
        virtual ~SLEPcHessian()
        {
            MatDestroy(&hessian);
            MatDestroy(&misforce);
            VecDestroy(&forcedipole_disp);
            EPSDestroy(&eps);
        };
        void setSystemData(std::shared_ptr<ParticleSystem> newsysdata)
        {
            m_sysdata = newsysdata;
        }
        void buildHessianandMisForce();
        
        void getMaxEigenvalue_forMumps();
        void checkSmallestEigenvalue_forMumps();
        void getAllEigenPairs_Mumps();
        void getEigenPairs();
        std::tuple<int,int> getPETScMatRange(unsigned int index); 
        
        double getEigenvalue(unsigned int index); 
        void saveEigenvalue(unsigned int index, std::shared_ptr<MPI::LogFile> logfile); 

        std::vector<double> getEigenvector(unsigned int index, bool forall); 
        void saveEigenvector(unsigned int index, std::shared_ptr<MPI::LogFile> logfile); 

        void calculateNonAffineTensor(); 
        void saveNonAffineTensor(unsigned int i, unsigned int j, std::shared_ptr<MPI::LogFile> logfile); 
        
        void solveForceDipoleProblem(double force_threshold);
        void saveForceDipoleProblem(std::shared_ptr<MPI::LogFile> logfile);
};

void SLEPcHessian::saveForceDipoleProblem(std::shared_ptr<MPI::LogFile> logfile)
{
    m_manager->notice(5) << "Saving computation" << std::endl;
    int globsize;
    Vec global_xr;
    VecScatter ctx;
    VecGetSize(forcedipole_disp,&globsize);
    VecScatterCreateToAll(forcedipole_disp,&ctx,&global_xr);
    
    double *tempvecp = new double[globsize];
    VecGetArray(global_xr,&tempvecp);
    VecScatterBegin(ctx,forcedipole_disp,global_xr,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,forcedipole_disp,global_xr,INSERT_VALUES,SCATTER_FORWARD);
    std::vector<double> displacement(tempvecp,tempvecp+globsize);//(hessian.rows());
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
            outline << displacement[2*id_x] << " ";
            outline << displacement[2*id_x+1] << " ";
            //end the input
            outline << std::endl;
            id_x+=1;
            logfile->write_shared(outline.str());
        }
    }
    //Destroy all objects
    VecRestoreArray(global_xr,&tempvecp);
    delete [] tempvecp;
    VecScatterDestroy(&ctx);
    VecDestroy(&global_xr);
    m_comm->barrier();
}

void SLEPcHessian::solveForceDipoleProblem(double force_threshold)
{
    Vec forcing; 
    KSP ksp;
    MatCreateVecs(hessian,NULL,&forcing);
    PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(hessian, &Istart, &Iend);
    //choose a random particle
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> int_dist(0,m_sysdata->particles.size()-1); // distribution in range [1, 6]
    bool particlenotfound = true;
    double forcedipole;
    PetscInt id_j, id_i;
    //All processes will look for
    m_manager->notice(5) << "Looking for pairs of particles to perturb" << std::endl; 
    while(particlenotfound)
    {
        //std::cout << dist6(rng) << std::endl;
        //py::print(id_i,id_i+1,Iend);
        int i = int_dist(rng);
        for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                    abr::get<position>(m_sysdata->particles[i]), max_rcut); p_j != false; ++p_j)
        {
            id_j = abr::get<abr::id>(*p_j);
            id_i = abr::get<abr::id>(m_sysdata->particles[i]);
            
            if (id_i != id_j && 2*id_i < Iend && 2*id_i >= Istart && 2*id_i + 1 < Iend && 2*id_i+1 >= Istart)
            {
                //Compute a list of local obsercavles 
                //Next we loop through the j-th particles for the virial stress
                //set the distance between particle i and particle j
                Eigen::Vector3d rij(-p_j.dx()[0], -p_j.dx()[1], -p_j.dx()[2]);
                m_potential->rij = rij;
                
                //Don't forget to set diameters of the potential
                m_potential->di =  abr::get<diameter>(m_sysdata->particles[i]);
                m_potential->dj =  abr::get<diameter>(*p_j);
                
                //Make sure the particle is unique
                if (m_potential->getRcut() > rij.dot(rij) && abs(m_potential->getPairForce()) < abs(force_threshold))
                {
                    particlenotfound = false;
                }
            }
        }
    }
    m_comm->barrier();

    //Now determine which rank has the smallest force. Once determined, we will broadcast the values
    std::vector< double > listofforces(m_comm->getSizeGlobal());
    m_comm->all_gather_v(forcedipole, listofforces);
    int therank = std::distance(listofforces.begin(),std::max_element(listofforces.begin(), listofforces.end()));
    m_comm->bcast(id_i, therank);
    m_comm->bcast(id_j, therank);
    
    m_manager->notice(5) << "Assembling forcing vector" << std::endl; 
    //m_comm->bcast(force, therank);
    for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                abr::get<position>(m_sysdata->particles[(int)id_i]), max_rcut); p_j != false; ++p_j)
    {
        id_j = abr::get<abr::id>(*p_j);
        id_i = abr::get<abr::id>(m_sysdata->particles[(int)id_i]);
        
        if (id_i != id_j && 2*id_i < Iend && 2*id_i >= Istart && 2*id_i + 1 < Iend && 2*id_i+1 >= Istart)
        {
            //Compute a list of local obsercavles 
            //Next we loop through the j-th particles for the virial stress
            //set the distance between particle i and particle j
            Eigen::Vector3d rij(-p_j.dx()[0], -p_j.dx()[1], -p_j.dx()[2]);
            m_potential->rij = rij;
            
            //Don't forget to set diameters of the potential
            m_potential->di =  abr::get<diameter>(m_sysdata->particles[(int)id_i]);
            m_potential->dj =  abr::get<diameter>(*p_j);
            
            //Make sure the particle is unique
            if (m_potential->getRcut() > rij.dot(rij) && abs(m_potential->getPairForce()) < abs(force_threshold))
            {
                m_manager->notice(5) << "Perturb Particle " << id_i << " and " << id_j << std::endl;
                forcedipole = m_potential->getPairForce();
                VecSetValue(forcing,2*id_i,-forcedipole*rij[0], ADD_VALUES);
                VecSetValue(forcing,2*id_i+1,-forcedipole*rij[1], ADD_VALUES);
                VecSetValue(forcing,2*id_j,forcedipole*rij[0], ADD_VALUES);
                VecSetValue(forcing,2*id_j+1,forcedipole*rij[1], ADD_VALUES);
                particlenotfound = false;
            }
        }
    }
    VecAssemblyBegin(forcing);
    VecAssemblyEnd(forcing);
    //set the nullspace
    m_manager->notice(5) << "Setting the null space" << std::endl; 
    MatNullSpace constant;
    MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,NULL,&constant);
    MatSetNullSpace(hessian,constant);

    //set the values according to particleid
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPSetOperators(ksp,hessian,hessian);
    KSPSetFromOptions(ksp);
    m_manager->notice(5) << "Solve!" << std::endl;
    KSPSolve(ksp,forcing,forcedipole_disp);
    MatNullSpaceDestroy(&constant);
};

SLEPcHessian::SLEPcHessian(std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, std::shared_ptr<PETScManager> manager, std::shared_ptr<MPI::Communicator> comm)
    : m_sysdata(sysdata), m_potential(potential), m_manager(manager), m_comm(comm)
    {
        ierr = PetscOptionsInsertString(NULL,m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
        m_manager->printPetscNotice(5,"Constructing SLEPc Hessian Object \n");
        
        double maxforce_rcut = potential->scaled_rcut;
        double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                                std::end(abr::get<diameter>(m_sysdata->particles)) );
        maxforce_rcut *= maxdiameter;
        max_rcut = std::max(maxforce_rcut,maxdiameter); //whichever is maximum
        
        //Construct Hessian
        m_manager->printPetscNotice(5,"Assemble PETSc Sparse Matrix \n");
        hessian_length = (unsigned int)m_sysdata->simbox->dim*m_sysdata->particles.size();
        
        buildHessianandMisForce();

        nconv = 0;
        m_maxeigval = 0;
        nonaffinetensor = Eigen::MatrixXd::Zero((unsigned int)m_sysdata->simbox->dim*((unsigned int)m_sysdata->simbox->dim+1)/2,(unsigned int)m_sysdata->simbox->dim*((unsigned int)m_sysdata->simbox->dim+1)/2);                
        MatCreateVecs(hessian,NULL,&forcedipole_disp);
        //Construct for loop here:
    };

void SLEPcHessian::buildHessianandMisForce()
{
        //MatCreate(PETSC_COMM_WORLD,&hessian);
        MatCreate(PETSC_COMM_WORLD,&hessian);
        MatSetSizes(hessian,PETSC_DETERMINE,PETSC_DETERMINE,hessian_length,hessian_length);
        MatSetType(hessian,MATAIJ);
        MatSetUp(hessian);

        //Construct mismatch force vector
        MatCreate(PETSC_COMM_WORLD,&misforce);
        MatSetSizes(misforce,PETSC_DETERMINE,PETSC_DETERMINE,hessian_length,(unsigned int)m_sysdata->simbox->dim*((unsigned int)m_sysdata->simbox->dim+1)/2);
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
                    m_potential->rij = rij;
                    
                    //Don't forget to set diameters of the potential
                    m_potential->di =  abr::get<diameter>(*p_i);
                    m_potential->dj =  abr::get<diameter>(*p_j);
                    
                       
                    //Make sure the particle is unique
                    if (m_potential->getRcut() > rij.dot(rij)) 
                    {
                        double factor = m_potential->getBondStiffness()+m_potential->getPairForce();
                        Eigen::Matrix3d offdiag_ij = -factor*nij*nij.transpose()
                                                     +Eigen::Matrix3d::Identity()*m_potential->getPairForce();

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

double SLEPcHessian::getEigenvalue(unsigned int index) 
{
        PetscReal lambda_r;//Vec xr,global_xr;
        ierr = EPSGetEigenvalue(eps,index,&lambda_r,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        return lambda_r;
};

void SLEPcHessian::saveEigenvalue(unsigned int index, std::shared_ptr<MPI::LogFile> logfile)
{
    double eigenval = getEigenvalue(index);
    if(m_comm->isRoot())
    {
        m_manager->notice(10) << "Process " << m_comm->getRank() << "is outputting eigenvalue" << std::endl;
        std::string outline = std::to_string(eigenval) + " ";
        logfile->write_shared(outline);
    }
}

std::vector<double> SLEPcHessian::getEigenvector(unsigned int index, bool forall) 
{
        int globsize;
        Vec xr,global_xr;
        VecScatter ctx;
        MatCreateVecs(hessian,NULL,&xr);
        VecGetSize(xr,&globsize);
        if (forall)
            VecScatterCreateToAll(xr,&ctx,&global_xr);
        else 
            VecScatterCreateToZero(xr,&ctx,&global_xr);
        
        double *tempvecp = new double[globsize];
        VecGetArray(global_xr,&tempvecp);
        ierr = EPSGetEigenvector(eps,index,xr,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        VecScatterBegin(ctx,xr,global_xr,INSERT_VALUES,SCATTER_FORWARD);
        VecScatterEnd(ctx,xr,global_xr,INSERT_VALUES,SCATTER_FORWARD);
        std::vector<double> eigenvector(tempvecp,tempvecp+globsize);//(hessian.rows());
        
        //Destroy all objects
        VecRestoreArray(global_xr,&tempvecp);
        delete [] tempvecp;
        VecScatterDestroy(&ctx);
        VecDestroy(&global_xr);
        VecDestroy(&xr);
        return eigenvector;
};

void SLEPcHessian::saveEigenvector(unsigned int index, std::shared_ptr<MPI::LogFile> logfile)
{
    std::vector<double> eigenvector = getEigenvector(index,true);
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

std::tuple<int,int> SLEPcHessian::getPETScMatRange(unsigned int index) 
{
        Vec xr;
        MatCreateVecs(hessian,NULL,&xr);
        int Istart,Iend;
        VecGetOwnershipRange(xr,&Istart, &Iend);
        VecDestroy(&xr);
        return std::make_tuple(Istart,Iend);
};

void SLEPcHessian::calculateNonAffineTensor() 
{
    if (nconv > 0)
    {
        ierr = EPSGetEigenvalue(eps,nconv-1,&m_maxeigval,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        double cond = m_manager->pinv_tol;//*tol;
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

void SLEPcHessian::saveNonAffineTensor(unsigned int i, unsigned int j, std::shared_ptr<MPI::LogFile> logfile)
{
    if(m_comm->isRoot())
    {
        m_manager->notice(10) << "Process " << m_comm->getRank() << "is outputting nonaffinetensor" << std::endl;
        std::string outline = std::to_string(nonaffinetensor(i,j)) + " ";
        logfile->write_shared(outline);
    }
}

void SLEPcHessian::getAllEigenPairs_Mumps()
{
    getMaxEigenvalue_forMumps();

    //Set value MUMPS work-space
    ierr = PetscOptionsInsertString(NULL,"-mat_mumps_icntl_14 200");CHKERRABORT(PETSC_COMM_WORLD,ierr);
    /*
     Set solver parameters at runtime
    */
    ierr = EPSSetFromOptions(eps);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    
    EPSType type;
    ST             st;          /* spectral transformation context */
    KSP            ksp;
    PC             pc;
    PetscInt       maxit;
    PetscReal      inta,intb,tol;
    
    ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    
    /*
        Spectrum slicing requires Krylov-Schur
    */
    ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSGetType(eps,&type);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Solution method: "+std::string(type)+"\n");
    
    ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Stopping condition: tol="+to_string_sci(tol)+", maxit="+std::to_string(maxit)+"\n");
    /*
        Set interval for spectrum slicing
    */
    inta = m_manager->lowerbound_tol;//*tol;
    intb = m_maxeigval*(1+inta);//PETSC_MAX_REAL;
    ierr = EPSSetInterval(eps,inta,intb);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Search interval is: ["+to_string_sci(inta)+", "+ to_string_sci(intb)+"]\n");
    
    /*
         Set shift-and-invert with Cholesky; select MUMPS if available
    */
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
        ierr = EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE);CHKERRABORT(PETSC_COMM_WORLD,ierr);  /* enforce zero detection */
        ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        /*
        Add several MUMPS options (see ex43.c for a better way of setting them in program):
        '-mat_mumps_icntl_13 1': turn off ScaLAPACK for matrix inertia
        '-mat_mumps_icntl_24 1': detect null pivots in factorization (for the case that a shift is equal to an eigenvalue)
        '-mat_mumps_cntl_3 <tol>': a tolerance used for null pivot detection (must be larger than machine epsilon)

        Note: depending on the interval, it may be necessary also to increase the workspace:
        '-mat_mumps_icntl_14 <percentage>': increase workspace with a percentage (50, 100 or more)
        */
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
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
     Get number of converged approximate eigenpairs
    */
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

void SLEPcHessian::getMaxEigenvalue_forMumps()
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
    m_manager->printPetscNotice(5,"Maximum Eigenvalue: "+to_string_sci(m_maxeigval)+"\n");
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
    m_manager->printPetscNotice(5,"Smallest Eigenvalue: "+to_string_sci(vr)+"\n");
};


void SLEPcHessian::getEigenPairs()
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
    /*
    Optional: Get some information from the solver and display it
    */
    ierr = EPSGetIterationNumber(eps,&its);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Number of iterations of the method: "+std::to_string(its)+"\n");
    ierr = EPSGetType(eps,&type);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Solution method: "+std::string(type)+"\n");
    ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Number of requested eigenvalues: "+std::to_string(nev)+"\n");
    ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_manager->printPetscNotice(5,"Stopping condition: tol="+to_string_sci(tol)+", maxit="+std::to_string(maxit)+"\n");
    
    ierr = MatCreateVecs(hessian,NULL,&xr);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = MatCreateVecs(hessian,NULL,&xi);CHKERRABORT(PETSC_COMM_WORLD,ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
     Get number of converged approximate eigenpairs
    */
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

void export_SLEPcHessian(py::module& m)
{
    py::class_<SLEPcHessian, std::shared_ptr<SLEPcHessian> >(m,"SLEPcHessian")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< PETScManager>  , std::shared_ptr< MPI::Communicator > >())
    .def("setSystemData", &SLEPcHessian::setSystemData)
    .def("buildHessianandMisForce", &SLEPcHessian::buildHessianandMisForce)
    .def("solveForceDipoleProblem", &SLEPcHessian::solveForceDipoleProblem)
    .def("saveForceDipoleProblem", &SLEPcHessian::saveForceDipoleProblem)
    .def("getEigenPairs", &SLEPcHessian::getEigenPairs)
    .def("getAllEigenPairs", &SLEPcHessian::getAllEigenPairs_Mumps)
    .def("calculateNonAffineTensor", &SLEPcHessian::calculateNonAffineTensor)
    .def("saveNonAffineTensor", &SLEPcHessian::saveNonAffineTensor)
    .def("saveEigenvector", &SLEPcHessian::saveEigenvector)
    .def("getEigenvalue", &SLEPcHessian::getEigenvalue)
    .def("saveEigenvalue", &SLEPcHessian::saveEigenvalue)
    .def("getPETScMatRange", &SLEPcHessian::getPETScMatRange)
    .def_readwrite("nconv", &SLEPcHessian::nconv)
    ;
};
#endif //__SLEPCHESSIAN_H__
