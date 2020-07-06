#ifndef __SLEPC_NMA_H__
#define __SLEPC_NMA_H__

#include "PETScHessianBase.h"
#include <slepceps.h>

class PYBIND11_EXPORT SLEPcNMA
{
    protected:
        std::shared_ptr< PETScHessianBase > m_hessian; 
        
        //List of observables to be computed!
        std::map<std::string, std::shared_ptr< PETScVectorFieldBase > > m_vectorfields;
        std::map<std::string, std::shared_ptr< PETScGlobalPropertyBase > > m_observables;
        
        int nconv; 
        double m_maxeigval;
        double m_mineigval;
        
        //PETSc objects for all computation 
        PetscErrorCode ierr; 
    
    public:
        SLEPcNMA(std::shared_ptr< PETScHessianBase > hessian);
        virtual ~SLEPcNMA()
        {
        };
        void setHessian(std::shared_ptr< PETScHessianBase > hessian)
        {
            m_hessian = hessian;
        }
        virtual void addVectorField(const std::shared_ptr< PETScVectorFieldBase >& obs)
        {
            m_vectorfields.insert(std::pair<std::string, std::shared_ptr< PETScVectorFieldBase > >(obs->name, obs));
        }

        virtual void addGlobalProperty(const std::shared_ptr< PETScGlobalPropertyBase >& obs)
        {
            m_observables.insert(std::pair<std::string, std::shared_ptr< PETScGlobalPropertyBase > >(obs->name, obs));
        }
        
        void getAllEigenPairs(std::string package);
        void getEigenPairs(std::string package);
        //Another function for getting some eigenpairs
        void getMaxEigenvalue();

        void saveResults(EPS& eps)
        {
            if (m_observables.count("nconv") > 0)
            {
                m_observables["nconv"]->setValue((double)nconv);
            }
            //Save Solution
            if(nconv >0)
            {
                calculateNonAffineTensor(eps);

                m_hessian->m_manager->printPetscNotice(5,"Storing selected eigenpairs \n");
                
                std::string prefix = "eigenvector_";
                std::string prefix1 = "eigenvalue_";
                std::string prefix2 = "eigenrelerror_";
                std::stringstream filename;
                for(int i = 0; i < nconv; ++i)
                {
                    filename.str("");
                    filename << prefix << i;
                    if (m_vectorfields.count(filename.str()) > 0)
                    {
                        m_vectorfields[filename.str()]->createVector(m_hessian->hessian);
                        ierr = EPSGetEigenvector(eps,i,m_vectorfields[filename.str()]->vectorobs,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                    }

                    filename.str("");
                    filename << prefix1 << i;
                    if (m_observables.count(filename.str()) > 0)
                    {
                        PetscReal vr;
                        ierr = EPSGetEigenvalue(eps,i,&vr,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                        m_observables[filename.str()]->setValue(vr);
                    }
                    
                    filename.str("");
                    filename << prefix2 << i;
                    if (m_observables.count(filename.str()) > 0)
                    {
                        PetscReal error;
                        ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error); CHKERRABORT(PETSC_COMM_WORLD,ierr);
                        m_observables[filename.str()]->setValue(error);
                    }
                }
            }
        }
        
        //void getEigenPairs();
        void calculateNonAffineTensor(const EPS& eps); 
        
        //double getEigenvalue(unsigned int index); 
        //std::vector<double> getEigenvector(unsigned int index, bool forall); 

        //void saveNonAffineTensor(unsigned int i, unsigned int j, std::shared_ptr<MPI::LogFile> logfile); 
};

SLEPcNMA::SLEPcNMA( std::shared_ptr< PETScHessianBase > hessian) 
    : m_hessian(hessian), nconv(0), m_maxeigval(std::numeric_limits<double>::max()), m_mineigval(-std::numeric_limits<double>::max())
{
    ierr = PetscOptionsInsertString(NULL,m_hessian->m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
};

void SLEPcNMA::getMaxEigenvalue()
{
    EPS eps;
    ST st;
    KSP ksp;
    PC pc;

    ierr = EPSCreate(PETSC_COMM_WORLD,&eps);
    ierr = EPSSetProblemType(eps,EPS_HEP); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSSetOperators(eps,m_hessian->hessian,NULL);
    
    m_hessian->m_manager->printPetscNotice(6,"Computing Maximum Eigenvalue \n");
    EPSSetDimensions(eps,1,PETSC_DEFAULT, PETSC_DEFAULT);
    EPSSetTolerances(eps,PETSC_DEFAULT,PETSC_DEFAULT); 
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);
    
    ierr = EPSGetST(eps,&st);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = STGetKSP(st,&ksp);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE);CHKERRABORT(PETSC_COMM_WORLD,ierr);  // enforce zero detection 
    ierr = PCFactorSetMatSolverType(pc,MATSOLVERPETSC);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    
    EPSSetUp(eps);
    EPSSolve(eps);
    
    ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSGetEigenvalue(eps,0,&m_maxeigval,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_hessian->m_manager->printPetscNotice(6,"Maximum Eigenvalue: "+detail::to_string_sci(m_maxeigval)+"\n");
    
    EPSDestroy(&eps);
};

void SLEPcNMA::getEigenPairs(std::string package)
{
    // Set solver parameters at runtime
    EPS eps;
    
    ierr = EPSCreate(PETSC_COMM_WORLD,&eps);
    ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSSetOperators(eps,m_hessian->hessian,NULL);
    ierr = EPSSetFromOptions(eps);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    if (package == "slepc-petsc")
    {
        ierr = PetscOptionsInsertString(NULL,"-st_pc_factor_mat_solver_type petsc");CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
    else if (package == "slepc-mumps")
    {
        #if defined(PETSC_HAVE_MUMPS)
            ierr = PetscOptionsInsertString(NULL,"-st_pc_factor_mat_solver_type mumps");CHKERRABORT(PETSC_COMM_WORLD,ierr);
        #endif
    }
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Solve the eigensystem
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ierr = EPSSetUp(eps);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_hessian->m_manager->printPetscNotice(5,"Solving for Eigenpairs . . . \n");
    ierr = EPSSolve(eps); CHKERRABORT(PETSC_COMM_WORLD,ierr);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                Display solution and clean up
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    //Get number of converged approximate eigenpairs
    ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_hessian->m_manager->printPetscNotice(5,"Number of converged eigenpairs: "+std::to_string(nconv)+"\n");
    if (m_hessian->m_manager->getNoticeLevel() >= 6)
    {
        m_hessian->m_manager->printPetscNotice(6,"Summary of Normal Mode Analysis: \n");
        ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
    
    saveResults(eps);
    EPSDestroy(&eps);
};

void SLEPcNMA::getAllEigenPairs(std::string package)
{
    if (m_hessian->areDiagonalsNonZero())
    {
        //Set the maximum eigenvalue
        getMaxEigenvalue();
        
        // Set solver parameters at runtime
        EPS eps;
        
        //A bunch of objects for solving eigenpairs
        EPSType type;
        ST             st;           
        KSP            ksp;
        PC             pc;
        PetscInt       maxit;
        PetscReal      inta,intb,tol;
       
        ierr = EPSCreate(PETSC_COMM_WORLD,&eps);
        ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSSetOperators(eps,m_hessian->hessian,NULL);
        ierr = EPSSetFromOptions(eps);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        
        //    Spectrum slicing requires Krylov-Schur
        ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSGetType(eps,&type);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        m_hessian->m_manager->printPetscNotice(5,"Solution method: "+std::string(type)+"\n");
        
        ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        m_hessian->m_manager->printPetscNotice(5,"Stopping condition: tol="+detail::to_string_sci(tol)+", maxit="+std::to_string(maxit)+"\n");
        //
        //   Set interval for spectrum slicing
        //
        inta = m_hessian->m_manager->lowerbound_tol; //tol;
        intb = m_maxeigval*(1+inta);//PETSC_MAX_REAL;
        ierr = EPSSetInterval(eps,inta,intb);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        m_hessian->m_manager->printPetscNotice(5,"Search interval is: ["+detail::to_string_sci(inta)+", "+ detail::to_string_sci(intb)+"]\n");
        
        //
        //     Set shift-and-invert with Cholesky; select PETSC default solver for this
        //
        ierr = EPSGetST(eps,&st);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = STSetType(st,STSINVERT);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = STGetKSP(st,&ksp);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = KSPSetType(ksp,KSPPREONLY);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = KSPGetPC(ksp,&pc);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        
        ierr = EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE);CHKERRABORT(PETSC_COMM_WORLD,ierr);  // enforce zero detection 
        if (package == "slepc-petsc")
        {
            ierr = PCSetType(pc,PCCHOLESKY);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PCFactorSetMatSolverType(pc,MATSOLVERPETSC);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PCFactorSetZeroPivot(pc,(PetscReal)std::numeric_limits<double>::epsilon()*m_hessian->m_manager->pivot_tol);
            ierr = PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PCFactorSetShiftAmount(pc, PETSC_DECIDE);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        }
        else if (package == "slepc-mumps")
        {
            #if defined(PETSC_HAVE_MUMPS)
                #if defined(PETSC_USE_COMPLEX)
                    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Spectrum slicing with MUMPS is not available for complex scalars");
                #endif
                ierr = PCSetType(pc,PCCHOLESKY);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                
                //Add several MUMPS options (see ex43.c for a better way of setting them in program):
                //'-mat_mumps_icntl_13 1': turn off ScaLAPACK for matrix inertia
                //'-mat_mumps_icntl_24 1': detect null pivots in factorization (for the case that a shift is equal to an eigenvalue)
                //'-mat_mumps_cntl_3 <tol>': a tolerance used for null pivot detection (must be larger than machine epsilon)

                //Note: depending on the interval, it may be necessary also to increase the workspace:
                //'-mat_mumps_icntl_14 <percentage>': increase workspace with a percentage (50, 100 or more)
                
                ierr = PetscOptionsInsertString(NULL,"-mat_mumps_icntl_13 1");CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = PetscOptionsInsertString(NULL,"-mat_mumps_icntl_24 1");CHKERRABORT(PETSC_COMM_WORLD,ierr);
                std::string cntl_3 = "-mat_mumps_cntl_3 ";
                cntl_3 +=  std::to_string(std::numeric_limits<float>::epsilon()*inta);
                ierr = PetscOptionsInsertString(NULL,cntl_3.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
            #endif
        }
        
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //                  Solve the eigensystem
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        ierr = EPSSetUp(eps);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        m_hessian->m_manager->printPetscNotice(5,"Solving for Eigenpairs . . . \n");
        ierr = EPSSolve(eps); CHKERRABORT(PETSC_COMM_WORLD,ierr);

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //                Display solution and clean up
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        //Get number of converged approximate eigenpairs
        ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        m_hessian->m_manager->printPetscNotice(5,"Number of converged eigenpairs: "+std::to_string(nconv)+"\n");
        if (m_hessian->m_manager->getNoticeLevel() >= 6)
        {
            m_hessian->m_manager->printPetscNotice(6,"Summary of Normal Mode Analysis: \n");
            ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
        }
        
        saveResults(eps);
        EPSDestroy(&eps);
    }
    else
    {
        m_hessian->m_manager->printPetscNotice(0,"[WARNING] Found particles with no neighbors during hessian assembly. Any normal mode analysis will be immediately aborted \n");
    }
};

void SLEPcNMA::calculateNonAffineTensor(const EPS& eps) 
{
    if (nconv > 0 && m_observables.count("nonaffinetensor") && m_hessian->areDiagonalsNonZero())
    {
        m_observables["nonaffinetensor"]->clear();
        m_hessian->m_manager->printPetscNotice(5,"Computing the nonaffine elasticity tensor \n");
        double cond = m_hessian->m_manager->pinv_tol;//tol;
        Vec xr, mult, seqmult;
        VecScatter ctx;
        
        PetscReal kr,re;
        MatCreateVecs(m_hessian->hessian,NULL,&xr);
        MatCreateVecs(m_hessian->misforce,&mult,NULL);
        
        //Create A Scatterrer and the resulting sequential vector holding temporary values 
        VecScatterCreateToAll(mult,&ctx,&seqmult);
        PetscInt colsize; MatGetSize(m_hessian->misforce,NULL,&colsize); 
        double *tempvecp = new double[colsize];//tempvec.data();
        VecGetArray(seqmult,&tempvecp);

        for (int I = 0; I < nconv; ++I)
        {   
            ierr = EPSGetEigenpair(eps,I,&kr,NULL,xr,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            #if defined(PETSC_USE_COMPLEX)
                re = PetscRealPart(kr);
            #else
                re = kr;
            #endif
            if (re > cond)
            {
                double laminv = 1/re;
                ierr = MatMultTranspose(m_hessian->misforce,xr,mult); CHKERRABORT(PETSC_COMM_WORLD,ierr);
                VecScatterBegin(ctx,mult,seqmult,INSERT_VALUES,SCATTER_FORWARD);
                VecScatterEnd(ctx,mult,seqmult,INSERT_VALUES,SCATTER_FORWARD);
                
                int Dim = m_hessian->m_sysdata->simbox->dim; 
                if (Dim == 2)
                {
                    for(int l = 0; l < m_observables["nonaffinetensor"]->getDimension(3); ++l)
                    {
                        for(int k = 0; k < m_observables["nonaffinetensor"]->getDimension(2); ++k)
                        {
                            for(int j = 0; j < m_observables["nonaffinetensor"]->getDimension(1); ++j)
                            {
                                for (int i = 0; i < m_observables["nonaffinetensor"]->getDimension(0); ++i)
                                {
                                    //if i,j or k,l  is 0,0 ---> 0
                                    //if i,j or k,l  is 0,1 ---> 1
                                    //if i,j or k,l  is 1,1 ---> 1
                                    int index = l+Dim*(k+Dim*(j+Dim*i));
                                    int id_i = i+j;
                                    int id_j = k+l;
                                    m_observables["nonaffinetensor"]->addValue(tempvecp[id_i]*tempvecp[id_j]*laminv,index);
                                }
                            }
                        }
                    }
                }
                m_hessian->m_comm->barrier();
            }
        }
        
        m_hessian->m_manager->printPetscNotice(5,"Deconstruct objects \n");
        VecRestoreArray(seqmult,&tempvecp);
        delete [] tempvecp;
        VecScatterDestroy(&ctx);
        VecDestroy(&seqmult);
        VecDestroy(&mult);
        VecDestroy(&xr);
        m_observables["nonaffinetensor"]->multiplyBy(1/m_hessian->m_sysdata->simbox->vol);
    }
    else
    {
        m_hessian->m_manager->printPetscNotice(0,"[WARNING] Found particles with no neighbors during hessian assembly. Any normal mode analysis will be immediately aborted \n");
    }
};

/*
double SLEPcNMA::getEigenvalue(unsigned int index) 
{
    PetscReal lambda_r;//Vec xr,global_xr;
    ierr = EPSGetEigenvalue(eps,index,&lambda_r,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    return lambda_r;
};

std::vector<double> SLEPcNMA::getEigenvector(unsigned int index, bool forall) 
{
    int globsize;
    Vec xr,global_xr;
    VecScatter ctx;
    MatCreateVecs(m_hessian->hessian,NULL,&xr);
    VecGetSize(xr,&globsize);
    std::vector<double> eigenvector;
    
    if (forall)
        VecScatterCreateToAll(xr,&ctx,&global_xr);
    else 
        VecScatterCreateToZero(xr,&ctx,&global_xr);
    
    double *tempvecp = new double[globsize];
    VecGetArray(global_xr,&tempvecp);
    VecScatterBegin(ctx,xr,global_xr,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx,xr,global_xr,INSERT_VALUES,SCATTER_FORWARD);
    if (forall)
    {
        std::copy(tempvecp, tempvecp+globsize, std::back_inserter(eigenvector));
    } 
    else
    {
        if (m_hessian->m_comm->isRoot())
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
*/

/*
template<int Dim>
void SLEPcNMA<Dim>::getEigenPairs()
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
        m_hessian->m_manager->printPetscNotice(5,"Performing spectrum slicing! Make sure all options are configured correctly \n");
    }
    
    EPSSetUp(eps);
    EPSSolve(eps);
    //Optional: Get some information from the solver and display it
    ierr = EPSGetIterationNumber(eps,&its);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_hessian->m_manager->printPetscNotice(5,"Number of iterations of the method: "+std::to_string(its)+"\n");
    ierr = EPSGetType(eps,&type);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_hessian->m_manager->printPetscNotice(5,"Solution method: "+std::string(type)+"\n");
    ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_hessian->m_manager->printPetscNotice(5,"Number of requested eigenvalues: "+std::to_string(nev)+"\n");
    ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_hessian->m_manager->printPetscNotice(5,"Stopping condition: tol="+detail::to_string_sci(tol)+", maxit="+std::to_string(maxit)+"\n");
    
    ierr = MatCreateVecs(m_hessian->hessian,NULL,&xr);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = MatCreateVecs(m_hessian->hessian,NULL,&xi);CHKERRABORT(PETSC_COMM_WORLD,ierr);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                Display solution and clean up
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    //
    // Get number of converged approximate eigenpairs
    //
    ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_hessian->m_manager->printPetscNotice(5,"Number of converged eigenpairs: "+std::to_string(nconv)+"\n");
    if (m_hessian->m_manager->getNoticeLevel() >= 6)
    {
        m_hessian->m_manager->printPetscNotice(6,"Summary of Normal Mode Analysis: \n");
        ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
    ierr = VecDestroy(&xr);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecDestroy(&xi);CHKERRABORT(PETSC_COMM_WORLD,ierr);
};

void SLEPcNMA::checkSmallestEigenvalue_forMumps()
{
    //This should only be done when computing a small amount of eigenvalues
    m_hessian->m_manager->printPetscNotice(5,"Computing Smallest Eigenvalue \n");
    //and setting the which to anything other than EPS_ALL
    EPSSetDimensions(eps,1,m_hessian->m_sysdata->particles.size(),m_hessian->m_sysdata->particles.size());
    EPSSetTolerances(eps,PETSC_DEFAULT,PETSC_DEFAULT); 
    EPSSetWhichEigenpairs(eps, EPS_SMALLEST_REAL);
    EPSSetUp(eps);
    EPSSolve(eps);
    
    PetscReal vr;
    ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSGetEigenvalue(eps,0,&vr,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_hessian->m_manager->printPetscNotice(5,"Smallest Eigenvalue: "+detail::to_string_sci(vr)+"\n");
    m_mineigval = vr;
    if (m_hessian->m_manager->getNoticeLevel() >= 6)
    {
        m_hessian->m_manager->printPetscNotice(6,"Summary of Normal Mode Analysis: \n");
        ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }
};
*/

void export_SLEPcNMA(py::module& m)
{
    py::class_<SLEPcNMA, std::shared_ptr<SLEPcNMA> >(m,"SLEPcNMA")
    .def(py::init< std::shared_ptr< PETScHessianBase > >())
    .def("setHessian", &SLEPcNMA::setHessian)
    .def("addVectorField", &SLEPcNMA::addVectorField)
    .def("addGlobalProperty", &SLEPcNMA::addGlobalProperty)
    .def("getAllEigenPairs", &SLEPcNMA::getAllEigenPairs)
    //.def("calculateNonAffineTensor", &T::calculateNonAffineTensor)
    //.def("getEigenPairs", &T::getEigenPairs)
    //.def("saveNonAffineTensor", &T::saveNonAffineTensor)
    //.def("saveEigenvector", &T::saveEigenvector)
    //.def("getEigenvalue", &T::getEigenvalue)
    //.def("saveEigenvalue", &T::saveEigenvalue)
    //.def("getPETScMatRange", &T::getPETScMatRange)
    ;
};

#endif //__SLEPC_NMA_H__
