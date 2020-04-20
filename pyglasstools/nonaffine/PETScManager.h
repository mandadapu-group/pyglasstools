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

class PYBIND11_EXPORT PETScManager : public Manager
{
    public:
        double pinv_tol;       
        PETScManager() : Manager(), ierr(0) 
        {
            auto argv = py::module::import("sys").attr("argv");
            pinv_tol = std::numeric_limits<double>::epsilon();
            for( auto test = argv.begin(); test != argv.end(); ++test)
            {
                if (std::string(py::str(*test)) == "-pinv_tol")
                {
                    //Increment manually
                    ++test;
                    try
                    {
                        pinv_tol = std::stod(py::str(*test),nullptr);
                    }
                    catch(...)
                    {
                        throw std::domain_error(std::string("Process [")+std::to_string(proc_rank)+std::string("]: ") +
                                                std::string("Did you correctly input your tolerance for pseudoinverse? The value you put is: ")+
                                                std::string(py::str(*test)));
                    }
                }
            }
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
    private:
        PetscErrorCode ierr; 
};

#endif
