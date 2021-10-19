#ifndef __PETSC_HESSIAN_BASE_H__
#define __PETSC_HESSIAN_BASE_H__

#include <pyglasstools/elasticity/HessianBase.h>
#include "PETScVectorField.h"
#include <petscksp.h>

/*
 * Derived base class for constructing Hessian using PETSc's matrices and arrays.
 * It can be further specialized for different applications, e.g., SLEPc and normal mode analysis.
 */
class PETScHessianBase : public HessianBase
{
    public:
        std::shared_ptr< PETScManager > m_manager; //!< a Manager class for printing values from PETSc arrays and outting error messages
        std::shared_ptr< MPI::ParallelCommunicator > m_comm;  //!< MPI Communicator object, wrapping common MPI methods
        Mat hessian; //!< the Hessian of the system, stored as a PETSc matrix
        Mat misforce; //!< the mismatch force vector, denoted as Xi in many papers. 
        PetscErrorCode ierr; //!< PETSC error code. A must-have for any calculations using PETSc. 

        PETScHessianBase(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                            std::shared_ptr< PETScManager > manager, std::shared_ptr< MPI::ParallelCommunicator > comm);
        virtual ~PETScHessianBase()
        {
        }
        
        virtual Mat getHessian()
        {
            return hessian;
        }
        
        virtual void setHessian(const Mat& inhessian)
        {
            hessian = inhessian;
        }
};

/*
 * Parametrized constructor. Args:
 * sysdata: the particle system, storing all of relevant configurational data
 * potential: the pair potential of the chosen system.
 * manager: a manager class for printing PETSc arrays and error messages.
 * comm: the MPI communicator
 *
 * All arguments are shared pointers, so that we are referring to an already-constructed objects. 
 */
PETScHessianBase::PETScHessianBase( std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                                    std::shared_ptr<PETScManager> manager, std::shared_ptr<MPI::ParallelCommunicator> comm)
    : HessianBase(sysdata,potential), m_manager(manager), m_comm(comm), ierr(0)
{
};

/*
 * Helper function to export PETScHessian base to Python
 */
void export_PETScHessianBase(py::module& m)
{
    py::class_<PETScHessianBase, HessianBase, std::shared_ptr<PETScHessianBase> >(m,"PETScHessianBase")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr<PETScManager>, std::shared_ptr< MPI::ParallelCommunicator > >())
    ;
};

#endif
