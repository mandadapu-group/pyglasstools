#ifndef __OBSERVABLES_H__
#define __OBSERVABLES_H__

#include <pyglasstools/MathAndTypes.h>
namespace abr = Aboria;

#include <pyglasstools/potential/PairPotential.h>

#include "../extern/pybind11/include/pybind11/pybind11.h"
namespace py = pybind11;
class PYBIND11_EXPORT Observable
{
    public:
        Observable() : name("NONE"), type("SCALAR"), islocal(true), isfield(false), useforce(false), dim(3) {};
        Observable(std::string _name, std::string _type, bool _islocal, bool _isfield, bool _useforce, int _dim) 
            : name(_name), type(_type), islocal(_islocal), isfield(_isfield), useforce(_useforce), dim(_dim)
            {
            };
        virtual void accumulate(const AboriaParticles::value_type& particle_i)
            {
                throw std::runtime_error("[ERROR] Observable does not require i-th particle data");
            }
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j)
            {
                throw std::runtime_error("[ERROR] Observable does not require i-th particle and/or j-th particle data");
            }
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j, 
                                const std::shared_ptr<PairPotential>& potential)
            {
                throw std::runtime_error("[ERROR] Observable does not require i-th particle, j-th particle, and/or potential force data");
            }
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                double cgval, unsigned int grid_id)
            {
                throw std::runtime_error("[error] observable is not a field");
            }
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j,
                                double bondval, unsigned int grid_id)
            {
                throw std::runtime_error("[error] observable is not a field");
            }
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j, 
                                const std::shared_ptr<PairPotential>& potential,
                                double bondval, unsigned int grid_id)
            {
                throw std::runtime_error("[error] observable is not a field");
            }
        virtual Eigen::MatrixXd getGlobalValue()
        {
            return Eigen::MatrixXd::Zero(1,1);
        }
        virtual std::vector< Eigen::MatrixXd > getField()
        {
            std::vector< Eigen::MatrixXd > val(1,Eigen::MatrixXd::Zero(1,1));
            return val;
        }
        std::string name;
        std::string type;
        bool islocal;
        bool isfield;
        bool useforce;
        int dim;
};


template< class AtomicObs >
class PYBIND11_EXPORT GlobalObservable : public Observable
{
    public:
        GlobalObservable() : Observable("NONE", "SCALAR", true, false, false,3), obs(3) {};
        GlobalObservable(std::string _name, std::string _type, bool _islocal, bool _useforce, int _dim) 
            : Observable(_name, _type, _islocal, false, _useforce, _dim), obs(_dim) 
            {
                if (type == "SCALAR")
                    val = Eigen::MatrixXd::Zero(1,1);
                else if (type == "VECTOR")
                    val = Eigen::MatrixXd::Zero(dim,1);
                else if (type == "TENSOR")
                    val = Eigen::MatrixXd::Zero(dim,dim);
                else
                    throw std::runtime_error("[ERROR] Type is unrecognized. Select from: SCALAR, VECTOR, and TENSOR");
            };
        virtual void accumulate(const AboriaParticles::value_type& particle_i)
            {
                obs.compute(particle_i, val);
            }
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j)
            {
                obs.compute(particle_i, particle_j, val);
            }
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j, 
                                const std::shared_ptr<PairPotential>& potential)
            {
                obs.compute(particle_i, particle_j, potential, val);
            }
        virtual Eigen::MatrixXd getGlobalValue()
        {
            return val;
        }
        
        AtomicObs obs;
        Eigen::MatrixXd val;
};

template< class AtomicObs >
class PYBIND11_EXPORT LocalObservable : public Observable
{
    public:
        LocalObservable() : Observable("NONE", "SCALAR", true, false, false,3), obs(3) {};
        LocalObservable(std::string _name, std::string _type, bool _islocal, bool _useforce, int _dim, int _gridsize) 
            : Observable(_name, _type, _islocal, false, _useforce, _dim), obs(_dim) 
            {
                val.resize(_gridsize);
                if (type == "SCALAR")
                    std::fill(val.begin(),val.end(),Eigen::MatrixXd::Zero(1,1));
                else if (type == "VECTOR")
                    std::fill(val.begin(),val.end(),Eigen::MatrixXd::Zero(dim,1));
                else if (type == "TENSOR")
                    std::fill(val.begin(),val.end(),Eigen::MatrixXd::Zero(dim,1));
                else
                    throw std::runtime_error("[ERROR] Type is unrecognized. Select from: SCALAR, VECTOR, and TENSOR");
            };
        virtual void accumulate(const AboriaParticles::value_type& particle_i,
                                double cgval, unsigned int grid_id)
            {
                obs.compute(particle_i, val[grid_id]);
                val[grid_id] *= cgval;
            }
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j,
                                double bondval, unsigned int grid_id)
            {
                obs.compute(particle_i, particle_j, val[grid_id]);
                val[grid_id] *= 0.5*bondval;
            }
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j, 
                                const std::shared_ptr<PairPotential>& potential,
                                double bondval, unsigned int grid_id)
            {
                obs.compute(particle_i, particle_j, potential, val[grid_id]);
                val[grid_id] *= 0.5*bondval;
            }

        virtual Eigen::MatrixXd getGlobalValue()
        {
            return std::accumulate(val.begin(), val.end(), Eigen::MatrixXd::Zero(dim,dim))/val.size();
        }
        virtual std::vector<Eigen::MatrixXd > getField()
        {
            return val;
        }
        
        AtomicObs obs;
        std::vector< Eigen::MatrixXd > val;
};

void export_Observable(py::module& m)
{
    py::class_<Observable, std::shared_ptr<Observable> >(m,"Observable")
    .def(py::init< std::string, std::string, bool, bool, bool, int >()) 
    //.def("accumulate", (void (T::*)(AboriaParticles::value_type)) &T::accumulate, "Accumulate value of local obs")
    //.def("accumulate", (void (T::*)(AboriaParticles::value_type,AboriaParticles::value_type)) &T::accumulate, "Accumulate value of pair obs")
    //.def("accumulate", (void (T::*)(AboriaParticles::value_type,AboriaParticles::value_type, const std::shared_ptr<PairPotential>&)) &T::accumulate, "Accumulate value of pair obs requiring force calc.")
    //.def("accumulate", (void (T::*)(AboriaParticles::value_type,double,unsigned int)) &T::accumulate, "Accumulate value of local obs")
    //.def("accumulate", (void (T::*)(AboriaParticles::value_type,AboriaParticles::value_type,std::shared_ptr<CoarseGrainFunction>)) &T::accumulate, "Accumulate value of pair obs")
    //.def("accumulate", (void (T::*)(AboriaParticles::value_type,AboriaParticles::value_type, const std::shared_ptr<PairPotential>&, std::shared_ptr<CoarseGrainFunction>)) &T::accumulate, "Accumulate value of pair obs requiring force calc.")
    .def("getGlobalValue",&Observable::getGlobalValue)
    .def("getField",&Observable::getField)
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
    .def_readwrite("val", &T::val)
    ;
};

/*
template<class T>
void export_LocalObservable(py::module& m, std::string name)
{
    py::class_<T, Observable, std::shared_ptr<T> >(m,name)
    .def(py::init< std::string, std::string, bool, bool, int >()) 
    .def_readwrite("val", &T::val)
    ;
};
*/
#endif
