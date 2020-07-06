#include "SimBox.h"

SimBox::SimBox( const double& _boxsize, 
                const Eigen::Vector3d& _origin, 
                const int& _dim)
{
    dim = _dim;
    origin = _origin;
    std::fill(boxsize.begin(),boxsize.end(),_boxsize);    
    setBounds();
    std::fill(periodic.begin(), periodic.end(),(int)1);    
    if (dim == 2)
    {
        boxsize[2] = 1.0;
        periodic[2] = 0;
        origin[2] = 0.0; //assert the zero in third axis of origin
        vol = boxsize[0]*boxsize[1];
    }
    else
    {
        vol = boxsize[0]*boxsize[1]*boxsize[2];
    }
};

SimBox::SimBox( const Eigen::Vector3d& _boxsize, 
                const Eigen::Vector3d& _origin, 
                const int& _dim)
{
    dim = _dim;
    origin = _origin;
    boxsize = _boxsize;//.fill(boxsize);    //.fill(0);
    setBounds();
    std::fill(periodic.begin(),periodic.end(),(int)1);    //.fill(0);
    if (dim == 2)
    {
        boxsize[2] = 1.0;
        periodic[2] = 0;
        origin[2] = 0.0; 
        vol = boxsize[0]*boxsize[1];
    }
    else
    {
        vol = boxsize[0]*boxsize[1]*boxsize[2];
    }
};

//! Sets the upperbound and lowerbound vector of the simulation box (what is the lowest vector contained in the box and vice versa
void SimBox::setBounds()
{
    upperbound = 0.5*boxsize-origin;
    lowerbound = (-1.0)*upperbound-2.0*origin;
}

//! Get the length of the box in each direction
void export_SimBox(py::module& m)
{
    py::class_<SimBox, std::shared_ptr<SimBox> >(m,"SimBox")
    .def(py::init<double, Eigen::Vector3d, int>())
    .def(py::init<Eigen::Vector3d, Eigen::Vector3d, int>())
    .def("minImage", &SimBox::minImage)   
    .def_readwrite("dim", &SimBox::dim)
    .def_readwrite("origin", &SimBox::origin)
    .def_readwrite("boxsize", &SimBox::boxsize)
    .def_readwrite("upperbound", &SimBox::upperbound)
    .def_readwrite("lowerbound", &SimBox::lowerbound)
    .def_readwrite("periodic", &SimBox::periodic)
    .def_readwrite("vol", &SimBox::vol)
    ;
};
