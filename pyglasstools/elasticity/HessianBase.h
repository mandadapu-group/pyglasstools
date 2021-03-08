#ifndef __PETSC_HESSIAN_BASE_H__
#define __PETSC_HESSIAN_BASE_H__

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
#include "GlobalProperty.h"
#include "PETScVectorField.h"

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
namespace py = pybind11;

#include <Aboria.h>
namespace abr = Aboria;
#include <petscksp.h>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/StdVector>

typedef Eigen::Triplet<double> Tripletd;
typedef Eigen::SparseMatrix<double> SparseMatd;

/*
 * Base class for constructing the Hessian of a system, given its pair potential and configuration data. 
 * Depending on which linear algebra package we use, then methods will be implemented differently. 
 */
class HessianBase
{
    public:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< PairPotential > m_potential; //!< system pari potential
        std::shared_ptr< MPI::Communicator > m_comm;  //!< MPI Communicator object, wrapping common MPI methods
        
        unsigned int hessian_length; //!< the total length of the Hessian matrix
        double max_rcut; //!<the maximum radius cut off for neighboring pairs of particles
        int diagonalsarenonzero; //!< check if all diagonal elements are non-zero!
       
        HessianBase(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                       std::shared_ptr< MPI::Communicator > comm);
        /* Default destructor */
        virtual ~HessianBase()
        {
        }
       
        /* 
         * Helper function to be used by the class destructor
         */ 
        virtual void destroyObjects()
        {
        };
        
        /* 
         * Helper function to be used by the class constructor
         */ 
        virtual void assembleObjects()
        {
        };
       
        /* 
         * Helper routine to check if the one of main diagonals of the Hessian matrix is non-zero or not.
         * This may happen if you have a particle which is trapped inside a cage but the interactions are
         * so short-ranged that it may-not feel any repulsive/attractive interactions from the surrounding particles.
         */
        bool areDiagonalsNonZero()
        {
            if (diagonalsarenonzero < 1)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
       
        /* 
         * Helper routine to set the system data. Useful when updating the currently stored data.
         */ 
        void setSystemData(std::shared_ptr<ParticleSystem> newsysdata)
        {
            m_sysdata = newsysdata;
        }
};

/*
 * Parametrized constructor. Args:
 * sysdata: the particle system, storing all of relevant configurational data
 * potential: the pair potential of the chosen system.
 * comm: the MPI communicator
 *
 * All arguments are shared pointers, so that we are referring to an already-constructed objects. 
 */
HessianBase::HessianBase( std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                          std::shared_ptr<MPI::Communicator> comm)
    : m_sysdata(sysdata), m_potential(potential), m_comm(comm), hessian_length(0), max_rcut(0.0)
{
};

/*
 * Helper function to export the base class HessianBase to Python
 */
void export_HessianBase(py::module& m)
{
    py::class_<HessianBase, std::shared_ptr<HessianBase> >(m,"HessianBase")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< MPI::Communicator > >())
    ;
};

/*
 * Derived base class for constructing Hessian using PETSc's matrices and arrays.
 * It can be further specialized for different applications, e.g., SLEPc and normal mode analysis.
 */
class PETScHessianBase : public HessianBase
{
    public:
        std::shared_ptr< PETScManager > m_manager; //!< a Manager class for printing values from PETSc arrays and outting error messages
        Mat hessian; //!< the Hessian of the system, stored as a PETSc matrix
        Mat misforce; //!< the mismatch force vector, denoted as Xi in many papers. 
        PetscErrorCode ierr; //!< PETSC error code. A must-have for any calculations using PETSc. 

        PETScHessianBase(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                            std::shared_ptr< PETScManager > manager, std::shared_ptr< MPI::Communicator > comm);
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
                                    std::shared_ptr<PETScManager> manager, std::shared_ptr<MPI::Communicator> comm)
    : HessianBase(sysdata,potential,comm), m_manager(manager), ierr(0)
{
};

/*
 * Helper function to export PETScHessian base to Python
 */
void export_PETScHessianBase(py::module& m)
{
    py::class_<PETScHessianBase, HessianBase, std::shared_ptr<PETScHessianBase> >(m,"PETScHessianBase")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr<PETScManager>, std::shared_ptr< MPI::Communicator > >())
    ;
};


/*
 * Derived base class for constructing Hessian using Eigen's matrices and arrays.
 * It can be further specialized for different applications, e.g., normal mode analysis.
 */
class EigenHessianBase : public HessianBase
{
    public:
        
        std::shared_ptr< EigenManager > m_manager; //!< a Manager class for printing values from Eigen arrays and outputing error messages
        SparseMatd hessian; //!< the Hessian of the system, stored as a sparse Eigen matrix
        Eigen::MatrixXd misforce; //!< the mismatch force vector, denoted as Xi in many papers. 

        EigenHessianBase(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                            std::shared_ptr< EigenManager > manager, std::shared_ptr< MPI::Communicator > comm);
        virtual ~EigenHessianBase()
        {
        }
        
};

/*
 * Parametrized constructor. Args:
 * sysdata: the particle system, storing all of relevant configurational data
 * potential: the pair potential of the chosen system.
 * manager: a manager class for printing Eigen arrays and error messages.
 * comm: the MPI communicator
 *
 * All arguments are shared pointers, so that we are referring to an already-constructed objects. 
 */
EigenHessianBase::EigenHessianBase( std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                                    std::shared_ptr<EigenManager> manager, std::shared_ptr<MPI::Communicator> comm)
    : HessianBase(sysdata,potential,comm), m_manager(manager)
{
};

/*
 * Helper function to export EigenHessian base to Python
 */
void export_EigenHessianBase(py::module& m)
{
    py::class_<EigenHessianBase, HessianBase, std::shared_ptr<EigenHessianBase> >(m,"EigenHessianBase")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr<EigenManager>, std::shared_ptr< MPI::Communicator > >())
    ;
};

#endif
