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
#include "PETScGlobalProperty.h"
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

class HessianBase
{
    public:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< PairPotential > m_potential; //!< system pari potential
        std::shared_ptr< MPI::Communicator > m_comm;  //!< MPI Communicator object, wrapping common MPI methods
        
        unsigned int hessian_length; //!< the total length of the Hessian matrix
        double max_rcut; //!<the maximum radius cut off for neighboring pairs of particles
        int diagonalsarenonzero; //!<check if all diagonal elements are non-zero!
        
        HessianBase(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                       std::shared_ptr< MPI::Communicator > comm);
        virtual ~HessianBase()
        {
        }
        
        virtual void destroyObjects()
        {
        };
        
        virtual void assembleObjects()
        {
        };
        
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
        
        void setSystemData(std::shared_ptr<ParticleSystem> newsysdata)
        {
            m_sysdata = newsysdata;
        }
};

HessianBase::HessianBase( std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                          std::shared_ptr<MPI::Communicator> comm)
    : m_sysdata(sysdata), m_potential(potential), m_comm(comm), hessian_length(0), max_rcut(0.0)
{
};

void export_HessianBase(py::module& m)
{
    py::class_<HessianBase, std::shared_ptr<HessianBase> >(m,"HessianBase")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< MPI::Communicator > >())
    ;
};

class PETScHessianBase : public HessianBase
{
    public:
        std::shared_ptr< PETScManager > m_manager;
        Mat hessian;
        Mat misforce;
        PetscErrorCode ierr;

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

PETScHessianBase::PETScHessianBase( std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                                    std::shared_ptr<PETScManager> manager, std::shared_ptr<MPI::Communicator> comm)
    : HessianBase(sysdata,potential,comm), m_manager(manager), ierr(0)
{
};

void export_PETScHessianBase(py::module& m)
{
    py::class_<PETScHessianBase, HessianBase, std::shared_ptr<PETScHessianBase> >(m,"PETScHessianBase")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr<PETScManager>, std::shared_ptr< MPI::Communicator > >())
    ;
};

class EigenHessianBase : public HessianBase
{
    public:
        std::shared_ptr< EigenManager > m_manager;
        SparseMatd hessian;
        Eigen::MatrixXd misforce;
        PetscErrorCode ierr;

        EigenHessianBase(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                            std::shared_ptr< EigenManager > manager, std::shared_ptr< MPI::Communicator > comm);
        virtual ~EigenHessianBase()
        {
        }
        
};

EigenHessianBase::EigenHessianBase( std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                                    std::shared_ptr<EigenManager> manager, std::shared_ptr<MPI::Communicator> comm)
    : HessianBase(sysdata,potential,comm), m_manager(manager), ierr(0)
{
};

void export_EigenHessianBase(py::module& m)
{
    py::class_<EigenHessianBase, HessianBase, std::shared_ptr<EigenHessianBase> >(m,"EigenHessianBase")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr<EigenManager>, std::shared_ptr< MPI::Communicator > >())
    ;
};

#endif
