//Look at SimBox.h
#ifndef __SIMBOX_H__
#define __SIMBOX_H__

#include <cmath>
#include "MathSTLVectors.h"

#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/stl.h"
namespace py = pybind11;

class PYBIND11_EXPORT SimBox
{
    public:
        //! Constructs a box from -m_boxsize/2 to m_boxsize/2
        /*! \param m_boxsize m_boxsizegth of one side of the box
            \post Box ranges from \c -m_boxsize/2 to \c +m_boxsize/2 in all 3 dimensions
            \post periodic = (1,1,1)
        */
        SimBox(const double& boxsize, const std::vector<double>& origin, const int& _dim)
        {
            dim = _dim;
            m_origin = origin;
            m_boxsize.resize(3);
            std::fill(m_boxsize.begin(),m_boxsize.end(),boxsize);    //.fill(0);
            setBounds();
            m_periodic.resize(3);
            std::fill(m_periodic.begin(),m_periodic.end(),(int)1);    //.fill(0);
            if (dim == 2)
            {
                m_boxsize[2] = 1.0;
                m_periodic[2] = 0;
                m_origin[2] = 0.0; //assert the zero in third axis of origin
            }
        };

        //! Constructs a box from -m_boxsize_x/2 to m_boxsize_x/2 for each dimension
        /*! \param m_boxsize_x m_boxsizegth of the x dimension of the box
            \param m_boxsize_y m_boxsizegth of the x dimension of the box
            \param m_boxsize_z m_boxsizegth of the x dimension of the box
            \post periodic = (1,1,1)
        */
        SimBox(const std::vector<double>& boxsize, const std::vector<double>& origin, const int& _dim)
        {
            dim = _dim;
            m_origin = origin;
            m_boxsize = boxsize;//.fill(boxsize);    //.fill(0);
            setBounds();
            m_periodic.resize(3);
            std::fill(m_periodic.begin(),m_periodic.end(),(int)1);    //.fill(0);
            if (dim == 2)
            {
                m_boxsize[2] = 1.0;
                m_periodic[2] = 0;
                m_origin[2] = 0.0; //assert the zero in third axis of origin
            }
        };
        
        ~SimBox(){};

        //! Set the periodic flags
        /*! \param periodic Flags to set
            \post Period flags are set to \a periodic
            \note It is invalid to set 1 for a periodic dimension where lo != -hi. This error is not checked for.
        */
        void setPeriodic(const std::vector<int>& periodic)
        {
            m_periodic = periodic;
        };
        
        //! Get the periodic flags
        /*! \return Periodic flags
        */
        const std::vector<int> getPeriodicVec()
        {
            return m_periodic;
        };
        bool getPeriodic(unsigned int i)
        {
            return (bool)m_periodic[i];
        };


        //! Update the box length
        /*! \param L new box length in each direction
        */
        void setBounds()
        {
            m_upperbound = 0.5*m_boxsize-m_origin;
            m_lowerbound = (-1.0)*m_upperbound-2.0*m_origin;
        }
        //! Get the length of the box in each direction
        /*! \returns The length of the box in each direction (hi - lo)
        */
        std::vector<double> getBoxSizeVec() const
        {
            return m_boxsize;
        }

        std::vector<double> getUpperBoundVec() const
        {
            return m_upperbound;
        }
        double getUpperBound(unsigned int i) const
        {
            return m_upperbound[i];
        }
        
        std::vector<double> getLowerBoundVec() const
        {
            return m_lowerbound;
        }
        double getLowerBound(unsigned int i) const
        {
            return m_lowerbound[i];
        }
        
        //! Update the box length
        /*! \param L new box length in each direction
        */
        void setDim(const int& _dim)
        {
            dim = _dim;
        }
        //! Get the length of the box in each direction
        /*! \returns The length of the box in each direction (hi - lo)
        */
        const int getDim()
        {
            return dim;
        }

        //! Get the volume of the box
        /*! \returns the volume
         *  \param twod If 1, return the area instead of the volume
         */
        const double getVolume()
        {
            if (dim == 2)
                return m_boxsize[0]*m_boxsize[1];
            else
                return m_boxsize[0]*m_boxsize[1]*m_boxsize[2];
        }

    protected:
        int dim;
        std::vector<double> m_origin;
        std::vector<double> m_boxsize;       //!< L precomputed (used to avoid subtractions in boundary conditions)
        std::vector<double> m_lowerbound;       //!< L precomputed (used to avoid subtractions in boundary conditions)
        std::vector<double> m_upperbound;       //!< minimum value of L, per coordinate precomputed
        std::vector<int> m_periodic; //!< 0/1 in each direction to tell if the box is periodic in that direction
        //Matrix<double, Dynamic, 3, ColMajor> m_gridposition;
};

// I might need to do something special since I'm using Eigen types
//an export function here
void export_SimBox(py::module& m)
{
    py::class_<SimBox, std::shared_ptr<SimBox> >(m,"SimBox")
    .def(py::init<double, std::vector<double>, int>())
    .def(py::init<std::vector<double>, std::vector<double>, int>())
    .def("getPeriodic", &SimBox::getPeriodic)
    .def("setPeriodic", &SimBox::setPeriodic)
    .def("getBoxSizeVec", &SimBox::getBoxSizeVec)
    .def("getUpperBoundVec", &SimBox::getUpperBoundVec)
    .def("getLowerBoundVec", &SimBox::getLowerBoundVec)
    .def("getDim", &SimBox::getDim)
    .def("setDim", &SimBox::setDim)
    .def("getVolume", &SimBox::getVolume)
    ;
};
#endif // __SIMBOX_H__
