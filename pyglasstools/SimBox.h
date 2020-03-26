//Look at SimBox.h
#ifndef __SIMBOX_H__
#define __SIMBOX_H__

#include <Eigen/Dense>
#include <cmath>
#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/eigen.h"
namespace py = pybind11;
using namespace Eigen;

class PYBIND11_EXPORT SimBox
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        //! Constructs a box from -Len/2 to Len/2
        /*! \param Len Length of one side of the box
            \post Box ranges from \c -Len/2 to \c +Len/2 in all 3 dimensions
            \post periodic = (1,1,1)
        */
        SimBox(double Len, int dim)
        {
            m_dim = dim;
            m_L = Vector3d::Constant(Len);    //.fill(0);
            setL();
            m_periodic = Vector3i::Constant(1);//fill(1);
            if (dim < 3)
                m_periodic[2] = 0;
        };

        //! Constructs a box from -Len_x/2 to Len_x/2 for each dimension
        /*! \param Len_x Length of the x dimension of the box
            \param Len_y Length of the x dimension of the box
            \param Len_z Length of the x dimension of the box
            \post periodic = (1,1,1)
        */
        SimBox(double Len_x, double Len_y, double Len_z, int dim)
        {
            m_dim = dim;
            m_L << Len_x, Len_y, Len_z;
            setL();
            m_periodic = Vector3i::Constant(1);//fill(1);
            if (dim < 3)
                m_periodic[2] = 0;
        };
        
        ~SimBox(){};

        //! Set the periodic flags
        /*! \param periodic Flags to set
            \post Period flags are set to \a periodic
            \note It is invalid to set 1 for a periodic dimension where lo != -hi. This error is not checked for.
        */
        void setPeriodic(const Vector3i& periodic)
        {
            m_periodic = periodic;
        };
        
        //! Get the periodic flags
        /*! \return Periodic flags
        */
        Vector3i getPeriodic() const
        {
            return m_periodic;
        };


        //! Update the box length
        /*! \param L new box length in each direction
        */
        void setL()
        {
            m_Lmax = m_L/2.0;//Scalar(2.0);
            m_Lmin = -m_Lmax;
        }
        //! Get the length of the box in each direction
        /*! \returns The length of the box in each direction (hi - lo)
        */
        Vector3d getL() const
        {
            return m_L;
        }

        Vector3d getLmax() const
        {
            return m_Lmax;
        }
        Vector3d getLmin() const
        {
            return m_Lmin;
        }
        
        //! Update the box length
        /*! \param L new box length in each direction
        */
        void setDim(const int &dim)
        {
            m_dim = dim;
        }
        //! Get the length of the box in each direction
        /*! \returns The length of the box in each direction (hi - lo)
        */
        int getDim() const
        {
            return m_dim;
        }

        //! Get the volume of the box
        /*! \returns the volume
         *  \param twod If 1, return the area instead of the volume
         */
        double getVolume()
        {
            if (m_dim == 2)
                return m_L[0]*m_L[1];
            else
                return m_L[0]*m_L[1]*m_L[2];
        }

    private:
        Matrix<double, Dynamic, 3, ColMajor> m_gridposition;
          
        Vector3d m_L;       //!< L precomputed (used to avoid subtractions in boundary conditions)
        Vector3d m_Lmin;       //!< L precomputed (used to avoid subtractions in boundary conditions)
        Vector3d m_Lmax;       //!< minimum value of L, per coordinate precomputed
        Vector3i m_periodic; //!< 0/1 in each direction to tell if the box is periodic in that direction
        int m_dim;
};

// I might need to do something special since I'm using Eigen types
//an export function here
void export_SimBox(py::module& m)
{
    py::class_<SimBox, std::shared_ptr<SimBox> >(m,"SimBox")
    .def(py::init<double, int>())
    .def(py::init<double, double, double, int>())
    .def("getPeriodic", &SimBox::getPeriodic)
    .def("setPeriodic", &SimBox::setPeriodic)
    .def("getL", &SimBox::getL)
    .def("setL", &SimBox::setL)
    .def("getDim", &SimBox::getDim)
    .def("setDim", &SimBox::setDim)
    .def("getVolume", &SimBox::getVolume)
    ;
};
#endif // __SIMBOX_H__
