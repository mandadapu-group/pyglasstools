#ifndef __OBSERVABLES_H__
#define __OBSERVABLES_H__
/*
#include <vector>
#include <cmath>
#include <string>
#include <map>

#include <algorithm>
#include <memory>
*/

#include <pyglasstools/MathAndTypes.h>
#include <pyglasstools/ParticleSystem.h>

#include "../extern/pybind11/include/pybind11/pybind11.h"
#include "../extern/pybind11/include/pybind11/stl.h"
#include "../extern/pybind11/include/pybind11/eigen.h"
namespace py = pybind11;

#include <Aboria.h>
namespace abr = Aboria;

class PYBIND11_EXPORT Observable
{
    public:
        Observable() : name("NONE"), type("SCALAR"), islocal(true), useforce(false), dim(3) {};
        Observable(std::string _name, std::string _type, bool _islocal, bool _useforce, int _dim) 
            : name(_name), type(_type), islocal(_islocal), useforce(_useforce), dim(_dim)
            {
                if (type == "SCALAR")
                    val = MatrixXd::Zero(1,1);//.resize(1,1);
                else if (type == "VECTOR")
                    val = MatrixXd::Zero(dim,1);//.resize(1,1);
                else if (type == "TENSOR")
                    val = MatrixXd::Zero(dim,dim);//.resize(1,1);
            };
        virtual void accumulate(AboriaParticles::value_type particle_i)//, AboriaParticles particle_j)
            {
            }
        virtual void accumulate(AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j)
            {
            }
        virtual void accumulate(AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j, 
                                const std::shared_ptr<PairPotential>& potential)
            {
            }
        
        std::string name;
        std::string type;
        bool islocal, useforce;
        int dim;
        MatrixXd val;
       
};

void export_Observable(py::module& m)
{
    py::class_<Observable, std::shared_ptr<Observable> >(m,"Observable")
    .def(py::init< std::string, std::string, bool, bool, int >()) 
    .def("accumulate", (void (Observable::*)(AboriaParticles::value_type)) &Observable::accumulate, "Accumulate value of local obs")
    .def("accumulate", (void (Observable::*)(AboriaParticles::value_type,AboriaParticles::value_type)) &Observable::accumulate, "Accumulate value of pair obs")
    .def("accumulate", (void (Observable::*)(AboriaParticles::value_type,AboriaParticles::value_type, const std::shared_ptr<PairPotential>&)) &Observable::accumulate, "Accumulate value of pair obs requiring force calc.")
    
    .def_readwrite("name", &Observable::name)
    .def_readwrite("type", &Observable::type)
    .def_readwrite("islocal", &Observable::islocal)
    .def_readwrite("useforce", &Observable::useforce)
    .def_readwrite("dim", &Observable::dim)
    .def_readwrite("val", &Observable::val)
    ;
};

#endif
