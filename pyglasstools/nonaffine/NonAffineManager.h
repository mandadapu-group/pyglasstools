#ifndef __PETSCMANAGER_H__
#define __PETSCMANAGER_H__

#include <pyglasstools/Manager.h>
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/embed.h>
namespace py = pybind11;

#include <petscsys.h>
#include <petscmat.h>
#include <petscvec.h>


class PYBIND11_EXPORT HessianManager : public Manager
{
    public:
        double pinv_tol, lowerbound_tol;  
        int nev, ncv, maxiter;
        std::string selrule;     
        
        HessianManager() : pinv_tol(1e-12), lowerbound_tol(1e-6), nev(1), ncv(2), maxiter(100), selrule("LM")
        {
            detail::argparser<double>("-pinv_tol",pinv_tol, "[ERROR] Invalid value for pinv_tol", argv);
            detail::argparser<double>("-lowerbound_tol",lowerbound_tol, "[ERROR] Invalid value for lowerbound_tol", argv);
            detail::argparser<int>("-spectra_nev",nev, "[ERROR] Invalid value spectra_nev", argv);
            detail::argparser<int>("-spectra_ncv",ncv, "[ERROR] Invalid value spectra_ncv", argv);
            detail::argparser<int>("-spectra_maxiter",maxiter, "[ERROR] Invalid value for maxier", argv);
            detail::argparser<std::string>("-spectra_selrule",selrule, "[ERROR] Invalid value for spectra selrule", argv);
        };
        ~HessianManager(){};
};


class PYBIND11_EXPORT PETScManager : public Manager
{
    public:
        PetscErrorCode ierr; 
        double pinv_tol, lowerbound_tol, upperbound_tol, pivot_tol, upperbound; 
        double fd_random_min, fd_random_max; 
        std::string fd_mode; 
        
        PETScManager() : ierr(0), pinv_tol(1e-12), lowerbound_tol(1e-6), upperbound_tol(1e-6), pivot_tol(1), upperbound(0.0), fd_random_min(0), fd_random_max(std::numeric_limits<double>::max()), fd_mode("random")  
        {
            detail::argparser<double>("-pinv_tol",pinv_tol, "[ERROR] Invalid value for pinv_tol", argv);
            detail::argparser<double>("-lowerbound_tol",lowerbound_tol, "[ERROR] Invalid value for lowerbound_tol", argv);
            detail::argparser<double>("-upperbound_tol",upperbound_tol, "[ERROR] Invalid value for upperbound_tol", argv);
            detail::argparser<double>("-upperbound",upperbound, "[ERROR] Invalid value for upperbound", argv);
            detail::argparser<double>("-pivot_tol", pivot_tol, "[ERROR] Invalid value for pivot_tol",argv);
            detail::argparser<double>("-fd_random_min", fd_random_min, "[ERROR] Invalid value for fd_random_min", argv);
            detail::argparser<double>("-fd_random_max", fd_random_max, "[ERROR] Invalid value for fd_random_max", argv);
            detail::argparser<std::string>("-fd_mode", fd_mode, "[ERROR] Invalid value for fd_mode",argv);
        };
        ~PETScManager(){};
        
        void printPetscNotice(unsigned int level,std::string statement)
        {
            if (level <= m_notice_level)
            {
                std::string prefix = m_notice_prefix+"("+std::to_string(level).c_str()+"): ";
                ierr = PetscPrintf(PETSC_COMM_WORLD,(prefix+statement).c_str());
                CHKERRABORT(PETSC_COMM_WORLD,ierr);
            }
        }           
        void printSynchronizedPetscNotice(unsigned int level, std::string statement)
        {
            if (level <= m_notice_level)
            {
                std::string prefix = m_notice_prefix+"("+std::to_string(level).c_str()+"): ";
                ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,(prefix+statement).c_str());
                CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
                CHKERRABORT(PETSC_COMM_WORLD,ierr);
            }
        }
        void viewMatrix(unsigned int level, const Mat& mat)
        {
            if (level <= m_notice_level)
            {
                std::string prefix = m_notice_prefix+"("+std::to_string(level).c_str()+"): ";
                ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,(prefix+std::string("View Matrix")).c_str());
                CHKERRABORT(PETSC_COMM_WORLD,ierr);
                ierr = MatView(mat, PETSC_VIEWER_STDOUT_WORLD);
                CHKERRABORT(PETSC_COMM_WORLD,ierr);
            }
        }
        void viewVector(unsigned int level, const Vec& vec)
        {
            if (level <= m_notice_level)
            {
                ierr = VecView(vec, PETSC_VIEWER_STDOUT_WORLD);
                CHKERRABORT(PETSC_COMM_WORLD,ierr);
            }
        }
};

void export_HessianManager(py::module& m)
{
    py::class_<HessianManager, std::shared_ptr<HessianManager> >(m,"HessianManager")
    .def(py::init<>())
    ;
};

void export_PETScManager(py::module& m)
{
    py::class_<PETScManager, std::shared_ptr<PETScManager> >(m,"PETScManager")
    .def(py::init<>())
    .def_readwrite("pinv_tol", &PETScManager::pinv_tol)
    .def_readwrite("lowerbound_tol", &PETScManager::lowerbound_tol)
    .def_readwrite("upperbound_tol", &PETScManager::upperbound_tol)
    ;
};

#endif
