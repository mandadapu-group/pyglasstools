#ifndef __OBSERVABLES_H__
#define __OBSERVABLES_H__

#include "MathAndTypes.h"
#ifdef ENABLE_MPI
#include "MPILogFile.h"
#else
#include "LogFile.h"
#endif
#include <pybind11/pybind11.h>
#include <pyglasstools/potential/PairPotential.h>

namespace py = pybind11;

//Public interface for storing observables
class PYBIND11_EXPORT Observable
{
    public:
        Observable() : name("NONE"), type("SCALAR"), islocal(true),  dim(3) {};
        Observable(std::string _name, std::string _type, bool _islocal, int _dim) 
            : name(_name), type(_type), islocal(_islocal), dim(_dim)
            {};
        virtual void clear(){};

        std::string name; //Name of the observable, e.g., pressure, etc.
        std::string type; //Type of the observable, e.g., the type is a vector, scalar, tensor, etc.
        bool islocal; //Does the observable require the property of a single particle, or does it require neighboring particles?
        int dim; //Dimensionality of the system
};


void export_Observable(py::module& m)
{
    py::class_<Observable, std::shared_ptr<Observable> >(m,"Observable")
    .def(py::init< std::string, std::string, bool, int >()) 
    .def_readwrite("name", &Observable::name)
    .def_readwrite("type", &Observable::type)
    .def_readwrite("islocal", &Observable::islocal)
    .def_readwrite("dim", &Observable::dim)
    ;
};

#endif
