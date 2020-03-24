//Look at SimBox.h
#ifndef __SIMBOX_H__
#define __SIMBOX_H__

#include <Eigen/Dense>
#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/eigen.h"
namespace py = pybind11;
using namespace Eigen;

class PYBIND11_EXPORT SimBox
{
    public:
        //A dummy constructor, just to have something setup when virtually no argument is given
        SimBox()
        {
            m_dim = 3; //default constructor will always give 3
            m_L = Vector3d::Zero(1.0);    //.fill(0);
            setL();
            m_periodic = Vector3i::Constant(1);//fill(1);
        };

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


        //! Apply periodic boundary conditions to a vector
        Vector3d applyPBC(const Vector3d v)
        {
            Vector3d w = v;
            //A bunch of if else statements here
            //z-direction
            if (m_periodic[2])
            {
                if (w[2] >= m_Lmax[2])
                    w[2] -= m_L[2];
                else if (w[2] < m_Lmin[2])
                    w[2] += m_L[2];
            }
            //y-direction
            if (m_periodic[1])
            {
                if (w[1] >= m_Lmax[1])
                    w[1] -= m_L[1];
                else if (w[1] < m_Lmin[1])
                    w[1] += m_L[1];
            }
            //x-direction
            if (m_periodic[0])
            {
                if (w[0] >= m_Lmax[0])
                {
                    w[0] -= m_L[0];
                    py::print(w[0]);
                }
                else if (w[0] < m_Lmin[0])
                {
                    w[0] += m_L[0];
                    py::print(w[0]);
                }
            }
            py::print(w[0],w[1],w[2]);
            return w;
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
    py::class_<SimBox>(m,"SimBox")
    .def(py::init<double, int>())
    .def(py::init<double, double, double, int>())
    .def("getPeriodic", &SimBox::getPeriodic)
    .def("setPeriodic", &SimBox::setPeriodic)
    .def("getL", &SimBox::getL)
    .def("setL", &SimBox::setL)
    .def("applyPBC", &SimBox::applyPBC)//minImage_overload)
    .def("getDim", &SimBox::getDim)
    .def("setDim", &SimBox::setDim)
    .def("getVolume", &SimBox::getVolume)
    ;
};
#endif // __SIMBOX_H__