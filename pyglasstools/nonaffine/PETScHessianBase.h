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

class PETScHessianBase
{
    public:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< PairPotential > m_potential; //!< system pari potential
        std::shared_ptr< PETScManager > m_manager; //!< manager object, which handlers error/logging/etc.
        std::shared_ptr< MPI::Communicator > m_comm;  //!< MPI Communicator object, wrapping common MPI methods
        
        Mat hessian;
        Mat misforce;
        unsigned int hessian_length;
        PetscErrorCode ierr;

        double max_rcut;
        
        PETScHessianBase(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                            std::shared_ptr< PETScManager > manager, std::shared_ptr< MPI::Communicator > comm);
        virtual ~PETScHessianBase()
        {
        }
        
        virtual void destroyPETScObjects()
        {
        };
        
        virtual void assemblePETScObjects()
        {
        };

        void setSystemData(std::shared_ptr<ParticleSystem> newsysdata)
        {
            m_sysdata = newsysdata;
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
                                    std::shared_ptr< PETScManager > manager, std::shared_ptr<MPI::Communicator> comm)
    : m_sysdata(sysdata), m_potential(potential), m_manager(manager), m_comm(comm)
{
};

void export_PETScHessianBase(py::module& m)
{
    py::class_<PETScHessianBase, std::shared_ptr<PETScHessianBase> >(m,"PETScHessianBase")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< PETScManager > , std::shared_ptr< MPI::Communicator > >())
    ;
};

#endif
