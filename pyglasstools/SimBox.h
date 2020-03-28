//Look at SimBox.h
#ifndef __SIMBOX_H__
#define __SIMBOX_H__

#include <cmath>
#include "MathAndTypes.h"

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
        SimBox(const double& _boxsize, const std::vector<double>& _origin, const int& _dim)
        {
            dim = _dim;
            origin = _origin;
            boxsize.resize(3);
            std::fill(boxsize.begin(),boxsize.end(),_boxsize);    //.fill(0);
            setBounds();
            periodic.resize(3);
            std::fill(periodic.begin(), periodic.end(),(int)1);    //.fill(0);
            if (dim == 2)
            {
                boxsize[2] = 1.0;
                periodic[2] = 0;
                origin[2] = 0.0; //assert the zero in third axis of origin
                vol = boxsize[0]*boxsize[1];
            }
            else
                vol = boxsize[0]*boxsize[1]*boxsize[2];
        };

        //! Constructs a box from -m_boxsize_x/2 to m_boxsize_x/2 for each dimension
        /*! \param m_boxsize_x m_boxsizegth of the x dimension of the box
            \param m_boxsize_y m_boxsizegth of the x dimension of the box
            \param m_boxsize_z m_boxsizegth of the x dimension of the box
            \post periodic = (1,1,1)
        */
        SimBox(const std::vector<double>& _boxsize, const std::vector<double>& _origin, const int& _dim)
        {
            dim = _dim;
            origin = _origin;
            boxsize = _boxsize;//.fill(boxsize);    //.fill(0);
            setBounds();
            periodic.resize(3);
            std::fill(periodic.begin(),periodic.end(),(int)1);    //.fill(0);
            if (dim == 2)
            {
                boxsize[2] = 1.0;
                periodic[2] = 0;
                origin[2] = 0.0; //assert the zero in third axis of origin
                vol = boxsize[0]*boxsize[1];
            }
            else
                vol = boxsize[0]*boxsize[1]*boxsize[2];
        };
        
        ~SimBox(){};

        //! Set the periodic flags
        /*! \param periodic Flags to set
            \post Period flags are set to \a periodic
            \note It is invalid to set 1 for a periodic dimension where lo != -hi. This error is not checked for.
        */
        void setPeriodic(const std::vector<int>& _periodic)
        {
            periodic = _periodic;
        };
        
        //! Get the periodic flags
        /*! \return Periodic flags
        */
        const std::vector<int> getPeriodicVec()
        {
            return periodic;
        };
        bool getPeriodic(unsigned int i)
        {
            return (bool)periodic[i];
        };


        //! Update the box length
        /*! \param L new box length in each direction
        */
        void setBounds()
        {
            upperbound = 0.5*boxsize-origin;
            lowerbound = (-1.0)*upperbound-2.0*origin;
        }
        //! Get the length of the box in each direction
        /*! \returns The length of the box in each direction (hi - lo)
        */
        std::vector<double> getBoxSizeVec() const
        {
            return boxsize;
        }

        std::vector<double> getUpperBoundVec() const
        {
            return upperbound;
        }
        double getUpperBound(unsigned int i) const
        {
            return upperbound[i];
        }
        
        std::vector<double> getLowerBoundVec() const
        {
            return lowerbound;
        }
        double getLowerBound(unsigned int i) const
        {
            return lowerbound[i];
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
                return boxsize[0]*boxsize[1];
            else
                return boxsize[0]*boxsize[1]*boxsize[2];
        }
        
        int dim;
        std::vector<double> origin;
        std::vector<double> boxsize;       //!< L precomputed (used to avoid subtractions in boundary conditions)
        std::vector<double> lowerbound;       //!< L precomputed (used to avoid subtractions in boundary conditions)
        std::vector<double> upperbound;       //!< minimum value of L, per coordinate precomputed
        std::vector<int> periodic; //!< 0/1 in each direction to tell if the box is periodic in that direction
        double vol;
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
    
    .def_readwrite("dim", &SimBox::dim)
    .def_readwrite("origin", &SimBox::origin)
    .def_readwrite("boxsize", &SimBox::boxsize)
    .def_readwrite("upperbound", &SimBox::upperbound)
    .def_readwrite("lowerbound", &SimBox::lowerbound)
    .def_readwrite("periodic", &SimBox::periodic)
    .def_readwrite("vol", &SimBox::vol)
    ;
};
#endif // __SIMBOX_H__
