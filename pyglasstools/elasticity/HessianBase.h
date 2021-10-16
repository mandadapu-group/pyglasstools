#ifndef __HESSIAN_BASE_H__
#define __HESSIAN_BASE_H__

#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>
#include <memory>

#include <pyglasstools/MathAndTypes.h>
#include <pyglasstools/ParticleSystem.h>
#include <pyglasstools/potential/PairPotential.h>

#include "GlobalProperty.h"

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
namespace py = pybind11;

#include <Aboria.h>
namespace abr = Aboria;

/*
 * Base class for constructing the Hessian of a system, given its pair potential and configuration data. 
 * Depending on which linear algebra package we use, then methods will be implemented differently. 
 */
class HessianBase
{
    public:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< PairPotential > m_potential; //!< system pari potential
        
        unsigned int hessian_length; //!< the total length of the Hessian matrix
        double max_rcut; //!<the maximum radius cut off for neighboring pairs of particles
        int diagonalsarenonzero; //!< check if all diagonal elements are non-zero!
       
        HessianBase(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential);
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
HessianBase::HessianBase( std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential)
    : m_sysdata(sysdata), m_potential(potential), hessian_length(0), max_rcut(0.0)
{
};

/*
 * Helper function to export the base class HessianBase to Python
 */
void export_HessianBase(py::module& m)
{
    py::class_<HessianBase, std::shared_ptr<HessianBase> >(m,"HessianBase")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential > >())
    ;
};

#endif
