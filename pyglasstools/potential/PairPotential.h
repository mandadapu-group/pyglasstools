#ifndef __PAIR_POTENTIAL_H__
#define __PAIR_POTENTIAL_H__

#include <Eigen/Dense>
#include <cmath>
#include "../extern/pybind11/include/pybind11/pybind11.h"
#include "../extern/pybind11/include/pybind11/eigen.h"
#include "../extern/pybind11/include/pybind11/stl.h"
namespace py = pybind11;
using namespace Eigen;

//Basic interface for common functions
class PYBIND11_EXPORT PairPotential
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        //Zero initialize everything
        PairPotential()
            : scaled_rcut(0), di(0), dj(0), rij(Vector3d::Zero())
        {
        };
        //If parameter is given
        PairPotential(const std::vector<double>& _params)
            : scaled_rcut(0), di(0), dj(0), rij(Vector3d::Zero()), params(_params)
        {
        };
        
        //If parameter and dimensionless cut-off radius is given
        PairPotential(const double& _scaled_rcut, const std::vector<double>& _params)
            : scaled_rcut(_scaled_rcut), di(0), dj(0), rij(Vector3d::Zero()), params(_params)
        {
        };
        virtual ~PairPotential(){};

        //Virtual methods here
        virtual double getRcut()
        {
            return 0.0;
        };
        //! Evaluate the force and energy
        virtual double getPairForce()
        {
            return 0.0;
        };
        
        double scaled_rcut;
        double di;
        double dj;
        Vector3d rij;
        std::vector<double> params;
};

template<class ShortRangeModel>
class PYBIND11_EXPORT ShortRangePairPotential : public PairPotential
{
    public:
        //Need this to avoid memory issue with STL Containers and Eigen Vectors
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        //Zero initialize everything
        ShortRangePairPotential(const std::vector<double>& _params)
            : PairPotential(_params)
        {
        };
        
        //Zero initialize everything
        ShortRangePairPotential(const double& _scaled_rcut, const std::vector<double>& _params)
            : PairPotential(_scaled_rcut, _params)
        {
        };
        ~ShortRangePairPotential(){};

        double getRcut()
        {
            double rsq_ij = rij.dot(rij);
            double rsq_cut = scaled_rcut*scaled_rcut;
            ShortRangeModel model(rsq_ij, rsq_cut, params);
            return model.computeRcut(di,dj);
        };
        //! Evaluate the force and energy
        double getPairForce()
        {
            double rsq_ij = rij.dot(rij);
            double rsq_cut = scaled_rcut*scaled_rcut;
            ShortRangeModel model(rsq_ij, rsq_cut, params);
            return model.computeForce(di,dj);
        };
};

void export_PairPotential(py::module& m)
{
    py::class_<PairPotential, std::shared_ptr<PairPotential> > (m, "PairPotential")
    .def(py::init<std::vector<double> >())
    .def(py::init<double, std::vector<double> >())
    .def("getRcut", &PairPotential::getRcut)
    .def("getPairForce", &PairPotential::getPairForce)
    .def_readwrite("scaled_rcut", &PairPotential::scaled_rcut)
    .def_readwrite("di", &PairPotential::di)
    .def_readwrite("dj", &PairPotential::dj)
    .def_readwrite("rij", &PairPotential::rij)
    .def_readwrite("params", &PairPotential::params)
    ;
};

template < class T > 
void export_ShortRangePairPotential(py::module& m, const std::string& name)
{
    py::class_<T, PairPotential, std::shared_ptr<T> >(m, name.c_str())
    .def(py::init<std::vector<double> >())
    .def(py::init<double, std::vector<double> >())
    ;
};




#endif // __PAIR_POTENTIAL_H__
