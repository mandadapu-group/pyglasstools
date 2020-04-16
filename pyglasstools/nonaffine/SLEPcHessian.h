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

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
namespace py = pybind11;

#include <Aboria.h>
namespace abr = Aboria;

#include <slepceps.h>
static char help[] = "Standard symmetric eigenproblem corresponding to the Laplacian operator in 2 dimensions.\n\n";

class PYBIND11_EXPORT SLEPcHessian
{
    public:
        Eigen::MatrixXd nonaffinetensor;
        int nconv; 
        double maxeigval;
        PetscReal upperbound;
        #pragma omp declare reduction( + : Eigen::MatrixXd : omp_out += omp_in ) initializer( omp_priv = Eigen::MatrixXd::Zero(omp_orig.rows(),omp_orig.cols()) )

        SLEPcHessian(std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential)
            : m_sysdata(sysdata), m_potential(potential)
            {
                auto argv = py::module::import("sys").attr("argv");
                int argc = 0;
                std::vector<const char*> cstrings;
                for( auto test = argv.begin(); test != argv.end(); ++test)
                {
                    std::string hey = py::str(*test);
                    cstrings.push_back(hey.data());
                    ++argc;
                }
                ierr = PetscPrintf(PETSC_COMM_WORLD,"Initialize SLEPc Object \n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
                SlepcInitialize(&argc,(char***)cstrings.data(),(char*)0,help);

                //py::module::import("sys").attr("argv");
                double maxforce_rcut = potential->scaled_rcut;
                double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                                        std::end(abr::get<diameter>(m_sysdata->particles)) );
                maxforce_rcut *= maxdiameter;
                max_rcut = std::max(maxforce_rcut,maxdiameter); //whichever is maximum
                
                //Construct Hessian
                ierr = PetscPrintf(PETSC_COMM_WORLD,"Assemble PETSc Sparse Matrix \n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
                hessian_length = (unsigned int)m_sysdata->simbox->dim*m_sysdata->particles.size();
                
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
                nconv = 0;
                maxeigval = 0;
                
                //MatView(misforce, PETSC_VIEWER_STDOUT_WORLD);
                //Construct for loop here:
                //Construct SLEPcHessian
            };
        virtual ~SLEPcHessian()
        {
            MatDestroy(&hessian);
            MatDestroy(&misforce);
            EPSDestroy(&eps);
            SlepcFinalize();
        };
        
        void calculateNonAffine(double tol) 
        {
            //nconv-1
            if (nconv > 0)
            {
                ierr = EPSGetEigenvalue(eps,nconv-1,&maxeigval,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                //double cond = hessian_length*maxeigval*tol;

                /*#pragma omp parallel
                { 
                    Vec xr,mult;
                    PetscReal kr,re;
                    MatCreateVecs(hessian,NULL,&xr);
                    MatCreateVecs(misforce,&mult,NULL);
                    
                    #pragma for reduction(+:nonaffinetensor) 
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
                            MatMultTranspose(misforce,xr,mult);
                        }
                    }
                }*/
                //nonaffinetensor /= m_sysdata->simbox->vol;
            }
        }

        void getAllEigenPairs_Mumps()
        {
            getMaxEigenvalue();
            //This should only be done when computing a small amount of eigenvalues
            //and setting the which to anything other than EPS_ALL
            EPSType type;
            ST             st;          /* spectral transformation context */
            KSP            ksp;
            PC             pc;
            PetscInt       nev,i,*inertias,ns,maxit;
            PetscReal      inta,intb,*shifts,tol;
            PetscBool      show=PETSC_FALSE;//,terse;
            
            //ierr = EPSSetTolerances(eps,_tol,_maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSSetWhichEigenpairs(eps,EPS_ALL);CHKERRABORT(PETSC_COMM_WORLD,ierr);

            /*
                Spectrum slicing requires Krylov-Schur
            */
            ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRABORT(PETSC_COMM_WORLD,ierr);

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
                ierr = PetscOptionsInsertString(NULL,"-mat_mumps_icntl_13 1 -mat_mumps_icntl_24 1 -mat_mumps_cntl_3 1e-50 -mat_mumps_icntl_14 100");CHKERRABORT(PETSC_COMM_WORLD,ierr);
            #endif

            /*
                Set interval for spectrum slicing
            */
            
            inta = 1e-8;
            intb = maxeigval+1;//PETSC_MAX_REAL;
            ierr = EPSSetInterval(eps,inta,intb);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            /*
             Set solver parameters at runtime
            */
            ierr = EPSSetFromOptions(eps);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSGetType(eps,&type);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Solution method: %s\n",type);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSGetInterval(eps,&inta,&intb);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Interval is  [%g, %g]\n",(double)inta,(double)intb);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //                  Solve the eigensystem
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ierr = EPSSetUp(eps);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            if (show) 
            {
                ierr = EPSKrylovSchurGetInertias(eps,&ns,&shifts,&inertias);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = PetscPrintf(PETSC_COMM_WORLD,"Subintervals (after setup):\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
                for (i=0;i<ns;i++) 
                { 
                    ierr = PetscPrintf(PETSC_COMM_WORLD,"Shift %g  Inertia %D \n",(double)shifts[i],inertias[i]);CHKERRABORT(PETSC_COMM_WORLD,ierr); 
                }
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = PetscFree(shifts);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = PetscFree(inertias);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            }
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Start Solving! \n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSSolve(eps); CHKERRABORT(PETSC_COMM_WORLD,ierr);
            
            ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            if (show) 
            {
            /*
             Show eigenvalues in interval and print solution
            */
            ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSGetInterval(eps,&inta,&intb);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"%D eigenvalues found in [%g, %g]\n",nev,(double)inta,(double)intb);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = EPSKrylovSchurGetInertias(eps,&ns,&shifts,&inertias);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = PetscPrintf(PETSC_COMM_WORLD,"All shifts (after solve):\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
                for (i=0;i<ns;i++) 
                { 
                    ierr = PetscPrintf(PETSC_COMM_WORLD,"Shift %g  Inertia %D \n",(double)shifts[i],inertias[i]);CHKERRABORT(PETSC_COMM_WORLD,ierr); 
                }
                ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = PetscFree(shifts);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = PetscFree(inertias);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            }
            ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
             
        };
        
        void getMaxEigenvalue()
        {
            //This should only be done when computing a small amount of eigenvalues
            //and setting the which to anything other than EPS_ALL
            EPSSetDimensions(eps,1,4,4);
            EPSSetTolerances(eps,PETSC_DEFAULT,PETSC_DEFAULT); 
            EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE);
            EPSSetUp(eps);
            EPSSolve(eps);
            
            ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of converged eigenpairs: %D\n",nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);

            ierr = EPSGetEigenvalue(eps,0,&maxeigval,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        };

        void getEigenPairs()
        {
            //This should only be done when computing a small amount of eigenvalues
            //and setting the which to anything other than EPS_ALL
            EPSType type;
            PetscReal      error,tol,re,im;
            PetscScalar    kr,ki;
            Vec            xr,xi;
            PetscInt       nev,maxit,its;
            
            EPSSetFromOptions(eps);
            EPSSetUp(eps);
            EPSWhich whicheig;
            EPSGetWhichEigenpairs(eps, &whicheig);
            if (whicheig == EPS_ALL)
            {
                ierr = PetscPrintf(PETSC_COMM_WORLD,"[WARNING] You're trying to compute eigenvalues in an interval (or the entirety of the spectrum)! \n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = PetscPrintf(PETSC_COMM_WORLD,"[WARNING] Please switch to eigsall instead of eigs function method! \n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
                //return;
            }
            EPSSolve(eps);
            /*
            Optional: Get some information from the solver and display it
            */
            ierr = EPSGetIterationNumber(eps,&its);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSGetType(eps,&type);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            
            ierr = MatCreateVecs(hessian,NULL,&xr);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = MatCreateVecs(hessian,NULL,&xi);CHKERRABORT(PETSC_COMM_WORLD,ierr);

            /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                            Display solution and clean up
             - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
            /*
             Get number of converged approximate eigenpairs
            */
            ierr = EPSGetConverged(eps,&nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = VecDestroy(&xr);CHKERRABORT(PETSC_COMM_WORLD,ierr);
            ierr = VecDestroy(&xi);CHKERRABORT(PETSC_COMM_WORLD,ierr);
        };

    protected:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< PairPotential > m_potential;
        double max_rcut;
        unsigned int hessian_length;
        
        Mat hessian;
        EPS eps;
        Mat misforce;
        PetscErrorCode ierr;
};

void export_SLEPcHessian(py::module& m)
{
    py::class_<SLEPcHessian, std::shared_ptr<SLEPcHessian> >(m,"SLEPcHessian")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >  >())
    .def("getEigenPairs", &SLEPcHessian::getEigenPairs)
    .def("getAllEigenPairs_Mumps", &SLEPcHessian::getAllEigenPairs_Mumps)
    //.def_readwrite("hessian", &SLEPcHessian::hessian)
    ;
};
#endif //__SLEPCHESSIAN_H__
