#ifndef __SLEPC_NMA_H__
#define __SLEPC_NMA_H__

#include "HessianBase.h"
#include <slepceps.h>


/*
 * A class for performing normal mode analysis and computing non-affine elasticity tensors
 * using the parallel eigensolver, SLEPc.
 * 
 * This class is the workhorse of the elasticity module. Any calculations of interest 
 * pertaining to normal mode analysis is contained within this class. 
 */
class PYBIND11_EXPORT SLEPcNMA
{
    protected:

        std::shared_ptr< PETScHessianBase > m_hessian; //<! pointer to the Hessian base class.
        
        std::map<std::string, std::shared_ptr< PETScVectorFieldBase > > m_vectorfields; //<! list of vector-field observables to compute
        std::map<std::string, std::shared_ptr< GlobalPropertyBase > > m_observables; //<! list of global observables to compute
        
        int nconv; //<! number of converged eigenmodes 
        double m_maxeigval; //<! maximum detected eigenvalue 
        double m_mineigval; //<! minimum detected eigenvalue
        
        //PETSc error code 
        PetscErrorCode ierr; 
    
    public:
        SLEPcNMA(std::shared_ptr< PETScHessianBase > hessian);
        virtual ~SLEPcNMA()
        {
        };
        
        /* Set the stored Hessian with another Hessian. Args:
         * hessian: a Hessian class 
         */
        void setHessian(std::shared_ptr< PETScHessianBase > hessian)
        {
            m_hessian = hessian;
        }

        /* Add a new vector field observable to a list of existing ones. Args:
         * obs: the new vector field observable
         */
        virtual void addVectorField(const std::shared_ptr< PETScVectorFieldBase >& obs)
        {
            m_vectorfields.insert(std::pair<std::string, std::shared_ptr< PETScVectorFieldBase > >(obs->name, obs));
        }
        
        /* Add a new global observable/property toa  list of existing ones. Args:
         * obs: the new global observable
         */
        virtual void addGlobalProperty(const std::shared_ptr< GlobalPropertyBase >& obs)
        {
            m_observables.insert(std::pair<std::string, std::shared_ptr< GlobalPropertyBase > >(obs->name, obs));
        }
       
        /* Obtain all possible eigenpairs, i.e., eigenvectors +eigenvalues. 
         * This is the routine that gets called whenever before we compute
         * the non-affine part of the elasticity tensor calculated. Args:
         * package: name of linear algebra package to use as a backend. 
         */ 
        void getAllEigenPairs(std::string package);

        /* Obtain subset of eigenpairs, i.e., eigenvectors +eigenvalues. Args:
         * package: name of linear algebra package to use as a backend.
         */ 
        void getEigenPairs(std::string package);
        
        /* Get the maximum eigenvalue of the system
         */
        void getMaxEigenvalue();
        
        /* Save the results of all computations. Gets called by 
         * getAllEigenPairs routine at the end. Args:
         * eps: class holding all computed eigenmodes
         */
        void saveResults(EPS& eps)
        {
            //To save the number of converged eigenmodes
            if (m_observables.count("nconv") > 0)
            {
                m_observables["nconv"]->setValue((double)nconv);
            }
            
            //To save and compute results dependent on eigenmodes
            if(nconv >0)
            {
                calculateNonAffineTensor(eps);

                m_hessian->m_manager->printPetscNotice(5,"Storing selected eigenpairs \n");
                
                //See if we want to save eigenvector, eigenvalue, and/or the relative error of the computation
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
        
        /* Computing the non-affine part of the elasticity tensor. Args:
         * eps: class holding all eigenmodes.
         */
        void calculateNonAffineTensor(const EPS& eps); 
};

/* Parametrized constructor. Args:
 * hessian: the input Hessian matrix
 */
SLEPcNMA::SLEPcNMA( std::shared_ptr< PETScHessianBase > hessian) 
    : m_hessian(hessian), nconv(0), m_maxeigval(std::numeric_limits<double>::max()), m_mineigval(-std::numeric_limits<double>::max())
{
    //Add any relevant command-line options for the PETSc linear algebra solver.
    ierr = PetscOptionsInsertString(NULL,m_hessian->m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
};


/* Get the maximum eigenvalue of the system
 */
void SLEPcNMA::getMaxEigenvalue()
{
    EPS eps; //<! the EPS object, holding all normal-mode analysis routines
    ST st; //<! the ST object, used for applying spectral transformation to the Hessian matrix
    KSP ksp; //<! the KSP object, holding basic linear algebra routines
    PC pc; //<! the PC object, holding informatino on pre-conditioners
    
    //Set up the EPS object
    ierr = EPSCreate(PETSC_COMM_WORLD,&eps);
    ierr = EPSSetProblemType(eps,EPS_HEP); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSSetOperators(eps,m_hessian->hessian,NULL);
    
    m_hessian->m_manager->printPetscNotice(6,"Computing Maximum Eigenvalue \n");

    //Maximum eigenvalue is efficiently computed using default setting of PETSc. 
    EPSSetDimensions(eps,1,PETSC_DEFAULT, PETSC_DEFAULT);
    EPSSetTolerances(eps,PETSC_DEFAULT,PETSC_DEFAULT); 
    EPSSetWhichEigenpairs(eps, EPS_LARGEST_REAL);
    
    //Additional setup for ST, KSP, and PC objects 
    ierr = EPSGetST(eps,&st);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = STGetKSP(st,&ksp);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE);CHKERRABORT(PETSC_COMM_WORLD,ierr);  // enforce zero detection 
    ierr = PCFactorSetMatSolverType(pc,MATSOLVERPETSC);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    
    //Final setup of the EPS object
    EPSSetUp(eps);
    
    //Solve for the maximum eigenvalue!
    EPSSolve(eps);
    
    ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = EPSGetEigenvalue(eps,0,&m_maxeigval,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    m_hessian->m_manager->printPetscNotice(6,"Maximum Eigenvalue: "+detail::to_string_sci(m_maxeigval)+"\n");
    
    EPSDestroy(&eps);
};

/* Obtain subset of eigenpairs, i.e., eigenvectors +eigenvalues. Args:
 * package: name of linear algebra package to use as a backend.
 */ 
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

/* Obtain all possible eigenpairs, i.e., eigenvectors +eigenvalues. 
 * This is the routine that gets called whenever before we compute
 * the non-affine part of the elasticity tensor calculated. Args:
 * package: name of linear algebra package to use as a backend. 
 */ 
void SLEPcNMA::getAllEigenPairs(std::string package)
{
    if (m_hessian->areDiagonalsNonZero())
    {
        //Set the maximum eigenvalue
        getMaxEigenvalue();
        
        // Set solver parameters at runtime
        EPS eps;
        
        //A bunch of PETSc objects for solving eigenpairs
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
        intb = m_maxeigval*(1+m_hessian->m_manager->upperbound_tol);//PETSC_MAX_REAL;
        if (m_hessian->m_manager->upperbound > 0.0)
        {
            intb = m_hessian->m_manager->upperbound*(1+m_hessian->m_manager->upperbound_tol);
        }
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
        
        //Depending on the solver backend, we will choose yet an additional set of settings to ensure numerics are done correctly
        if (package == "slepc-petsc")
        {
            ierr = PCSetType(pc,PCCHOLESKY);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PCFactorSetMatSolverType(pc,MATSOLVERPETSC);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PCFactorSetZeroPivot(pc,(PetscReal)std::numeric_limits<float>::epsilon()*m_hessian->m_manager->pivot_tol);
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
                cntl_3 +=  std::to_string(std::numeric_limits<float>::epsilon()*m_hessian->m_manager->pivot_tol);
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


/* Computing the non-affine part of the elasticity tensor. Args:
 * eps: class holding all eigenmodes.
 */
void SLEPcNMA::calculateNonAffineTensor(const EPS& eps) 
{
    if (nconv > 0 && m_observables.count("nonaffinetensor") && m_hessian->areDiagonalsNonZero())
    {
        m_hessian->m_manager->printPetscNotice(5,"Computing the nonaffine elasticity tensor \n");
        
        m_observables["nonaffinetensor"]->clear();
        double cond = m_hessian->m_manager->pinv_tol;
        Vec xr, mult, seqmult;
        VecScatter ctx;
        
        PetscReal kr,re;

        //We need a bit of prep because the layout may or may not be compatible between hessian and misforce
        MatCreateVecs(m_hessian->hessian,NULL,&xr);
        MatCreateVecs(m_hessian->misforce,&mult,NULL);
        
        //In addition, we need to a temporary array, tempvecp, to store the resulting dot product between mismatch force vector
        //and the eigenvector. 
        VecScatterCreateToAll(mult,&ctx,&seqmult);
        PetscInt colsize; MatGetSize(m_hessian->misforce,NULL,&colsize); 
        double *tempvecp = new double[colsize];
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
                //Dot product of mismatch force vector with an eigenvector
                ierr = MatMultTranspose(m_hessian->misforce,xr,mult); CHKERRABORT(PETSC_COMM_WORLD,ierr);
                
                //Scatter and gather the resulting dot products, so that each process actually possess all values
                VecScatterBegin(ctx,mult,seqmult,INSERT_VALUES,SCATTER_FORWARD);
                VecScatterEnd(ctx,mult,seqmult,INSERT_VALUES,SCATTER_FORWARD);
                
                int Dim = m_hessian->m_sysdata->simbox->dim; 
                
                //Although the non-affine tensor is stores as a proper 4-th rank tensor, the mismatch force vector, i.e., misforce
                //is stored as a flattened array. There is a specific mapping between the flattened array to the 4-th rank tensor. 
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
                                if (Dim == 3 && (id_i == 2 || id_j == 2))
                                {
                                    if ((i == 0 && j == 2) || (i == 2 && j == 0) )
                                    {
                                        id_i = 5;
                                    }
                                    if ((k == 0 && l == 2) || (k == 2 && l == 0))
                                    {
                                        id_j = 5;
                                    }
                                    m_observables["nonaffinetensor"]->addValue(tempvecp[id_i]*tempvecp[id_j]*laminv,index);
                                }
                                else
                                {
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

        //Don't forget to divide by the system volume!
        m_observables["nonaffinetensor"]->multiplyBy(1/m_hessian->m_sysdata->simbox->vol);
    }
    else
    {
        if (nconv <= 0)
        {
            m_hessian->m_manager->printPetscNotice(0,"[WARNING] No converged eigenpairs available \n");
        }
        if (m_observables.count("nonaffinetensor") == 0)
        {
            m_hessian->m_manager->printPetscNotice(0,"[WARNING] Not calculating non-affine elasticity tensor \n");
        }
        if (!m_hessian->areDiagonalsNonZero())
        {
            m_hessian->m_manager->printPetscNotice(0,"[WARNING] Found particles with no neighbors during hessian assembly. Any normal mode analysis will be immediately aborted \n");
        } 
    }
};

/*
 * Helper function to exprot SLEPcNMA to Python
 */
void export_SLEPcNMA(py::module& m)
{
    py::class_<SLEPcNMA, std::shared_ptr<SLEPcNMA> >(m,"SLEPcNMA")
    .def(py::init< std::shared_ptr< PETScHessianBase > >())
    .def("setHessian", &SLEPcNMA::setHessian)
    .def("addVectorField", &SLEPcNMA::addVectorField)
    .def("addGlobalProperty", &SLEPcNMA::addGlobalProperty)
    .def("getAllEigenPairs", &SLEPcNMA::getAllEigenPairs)
    ;
};

#endif //__SLEPC_NMA_H__
