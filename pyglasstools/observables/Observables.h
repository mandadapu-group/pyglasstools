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
//#include <pyglasstools/cgfunc/CoarseGrainFunction.h>

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
            };
        virtual void accumulate(AboriaParticles::value_type particle_i)
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
        virtual MatrixXd getGlobalValue()
        {
            return MatrixXd::Zero(1,1);
        }
        std::string name;
        std::string type;
        bool islocal;
        bool useforce;
        int dim;
};


template< class AtomicObs >
class PYBIND11_EXPORT GlobalObservable : public Observable
{
    public:
        GlobalObservable() : Observable("NONE", "SCALAR", true, false,3), obs(3) {};
        GlobalObservable(std::string _name, std::string _type, bool _islocal, bool _useforce, int _dim) 
            : Observable(_name, _type, _islocal, _useforce, _dim), obs(_dim) 
            {
                if (type == "SCALAR")
                    val = MatrixXd::Zero(1,1);
                else if (type == "VECTOR")
                    val = MatrixXd::Zero(dim,1);
                else if (type == "TENSOR")
                    val = MatrixXd::Zero(dim,dim);
                else
                    throw std::runtime_error("[ERROR] Type is unrecognized. Select from: SCALAR, VECTOR, and TENSOR");
            };
        virtual void accumulate(AboriaParticles::value_type particle_i)
            {
                obs.compute(particle_i, val);
            }
        virtual void accumulate(AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j)
            {
                obs.compute(particle_i, particle_j, val);
            }
        virtual void accumulate(AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j, 
                                const std::shared_ptr<PairPotential>& potential)
            {
                obs.compute(particle_i, particle_j, potential, val);
            }
        virtual MatrixXd getGlobalValue()
        {
            return val;
        }
        
        AtomicObs obs;
        MatrixXd val;
};

/*
template< class AtomicObs >
class PYBIND11_EXPORT LocalObservable : public Observable
{
    public:
        LocalObservable() : Observable("NONE", "SCALAR",true,false,3) {};
        LocalObservable(std::string _name, std::string _type, bool _islocal, bool _useforce, int _dim, int _gridsize) 
            : Observable(_name, _type, _islocal, _useforce, _dim) 
            {
                val.resize(_gridsize);
                if (type == "SCALAR")
                    std::fill(val.begin(),val.begin(),MatrixXd::Zero(1,1));
                else if (type == "VECTOR")
                    std::fill(val.begin(),val.begin(),MatrixXd::Zero(dim,1));
                else if (type == "TENSOR")
                    std::fill(val.begin(),val.begin(),MatrixXd::Zero(dim,dim));
                else
                    throw std::runtime_error("[ERROR] Type is unrecognized. Select from: SCALAR, VECTOR, and TENSOR");
            };
        virtual void accumulate(AboriaParticles::value_type particle_i,
                                std::shared_ptr<CoarseGrainFunction> cgfunc,
                                unsigned int grid_id)
            {
                AtomicObs obs(particle_i,cgfunc,grid_id)
                obs.compute(val);
            }
        virtual void accumulate(AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j,
                                std::shared_ptr<CoarseGrainFunction> cgfunc, 
                                unsigned int grid_id)
            {
                AtomicObs obs(particle_i,particle_j, cgfunc,grid_id)
                obs.compute(val);
            }
        virtual void accumulate(AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j, 
                                const std::shared_ptr<PairPotential>& potential,
                                std::shared_ptr<CoarseGrainFunction>, unsigned int grid_id)
            {
                AtomicObs obs(particle_i,particle_j, potential, cgfunc,grid_id)
                obs.compute(val);
            }

        std::vector< MatrixXd > val;
};
*/

void export_Observable(py::module& m)
{
    py::class_<Observable, std::shared_ptr<Observable> >(m,"Observable")
    .def(py::init< std::string, std::string, bool, bool, int >()) 
    .def("getGlobalValue",&Observable::getGlobalValue)
    .def_readwrite("name", &Observable::name)
    .def_readwrite("type", &Observable::type)
    .def_readwrite("islocal", &Observable::islocal)
    .def_readwrite("useforce", &Observable::useforce)
    .def_readwrite("dim", &Observable::dim)
    ;
};

template<class T>
void export_GlobalObservable(py::module& m, const std::string& name)
{
    py::class_<T, Observable, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::string, std::string, bool, bool, int >()) 
    .def("accumulate", (void (T::*)(AboriaParticles::value_type)) &T::accumulate, "Accumulate value of local obs")
    .def("accumulate", (void (T::*)(AboriaParticles::value_type,AboriaParticles::value_type)) &T::accumulate, "Accumulate value of pair obs")
    .def("accumulate", (void (T::*)(AboriaParticles::value_type,AboriaParticles::value_type, const std::shared_ptr<PairPotential>&)) &T::accumulate, "Accumulate value of pair obs requiring force calc.")
    .def_readwrite("val", &T::val)
    ;
};

/*
template<class T>
void export_LocalObservable(py::module& m, std::string name)
{
    py::class_<T, Observable, std::shared_ptr<T> >(m,name)
    .def(py::init< std::string, std::string, bool, bool, int >()) 
    .def("accumulate", (void (T::*)(AboriaParticles::value_type,std::shared_ptr<CoarseGrainFunction>)) &T::accumulate, "Accumulate value of local obs")
    .def("accumulate", (void (T::*)(AboriaParticles::value_type,AboriaParticles::value_type,std::shared_ptr<CoarseGrainFunction>)) &T::accumulate, "Accumulate value of pair obs")
    .def("accumulate", (void (T::*)(AboriaParticles::value_type,AboriaParticles::value_type, const std::shared_ptr<PairPotential>&, std::shared_ptr<CoarseGrainFunction>)) &T::accumulate, "Accumulate value of pair obs requiring force calc.")
    .def_readwrite("val", &T::val)
    ;
};
*/
#endif
