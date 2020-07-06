#ifndef __SIMBOX_H__
#define __SIMBOX_H__

#include <cmath>
#include "MathAndTypes.h"

#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

class PYBIND11_EXPORT SimBox
{
    public:
        int dim;

        Eigen::Vector3d origin;         //!< origin vector
        Eigen::Vector3d boxsize;       //!< L precomputed (used to avoid subtractions in boundary conditions)
        Eigen::Vector3d lowerbound;       //!< L precomputed (used to avoid subtractions in boundary conditions)
        Eigen::Vector3d upperbound;       //!< minimum value of L, per coordinate precomputed
        Eigen::Vector3d periodic; //!< 0/1 in each direction to tell if the box is periodic in that direction
        
        double vol;
        
        //! Constructs a box from -m_boxsize/2 to m_boxsize/2
        /*! \param _boxsize length of one side of simulation box. Avery side is assumed to have the same length
            \param _origin the origin of the box, which may differ from (0,0,0)
            \param _dim dimensionality of the box
        */
        SimBox( const double& _boxsize, 
                const Eigen::Vector3d& _origin, 
                const int& _dim);
        //! Constructs a box from -m_boxsize_x/2 to m_boxsize_x/2 for each dimension
        /*! \param _boxsize a vector of length for each side of simulation box.
            \param _origin the origin of the box, which may differ from (0,0,0)
            \param _dim dimensionality of the box
        */
        SimBox( const Eigen::Vector3d& _boxsize, 
                const Eigen::Vector3d& _origin, 
                const int& _dim);
       
        //! Empty destructor 
        virtual ~SimBox(){};

        //! Sets the upperbound and lowerbound vector of the simulation box (what is the lowest vector contained in the box and vice versa
        void setBounds();
        
        Eigen::Vector3d minImage(const Eigen::Vector3d& v) const
        {
            Eigen::Vector3d w = v;
            if (periodic[2])
            {
                double img = rintf(w[2] / boxsize[2]);
                w[2] -= boxsize[2] * img;
            }
            if (periodic[1])
            {
                double img = rintf(w[1] / boxsize[1]);
                w[1] -= boxsize[1] * img;
            }
            if (periodic[0])
            {
                double img = rintf(w[0] / boxsize[0]);
                w[0] -= boxsize[0] * img;
            }
            return w;
        }
};

// I might need to do something special since I'm using Eigen types
//an export function here
void export_SimBox(py::module& m);

#endif // __SIMBOX_H__
