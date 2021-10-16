#ifndef __PAIR_POTENTIAL_H__
#define __PAIR_POTENTIAL_H__

#include <Eigen/Dense>
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Eigen;

//Basic interface for common functions
class PYBIND11_EXPORT PairPotential
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        //Zero initialize everything
        PairPotential()
            : scaled_rcut(0)
        {
        };
        //If parameter is given
        PairPotential(const std::vector<double>& _params)
            : scaled_rcut(0), params(_params)
        {
        };
        
        //If parameter and dimensionless cut-off radius is given
        PairPotential(const double& _scaled_rcut, const std::vector<double>& _params)
            : scaled_rcut(_scaled_rcut), params(_params)
        {
        };
        virtual ~PairPotential(){};

        virtual double getRcut(Eigen::Vector3d rij, double di, double dj)
        {
           return 0;
        }
        
        virtual void setParams(double val, int num)
        {
            params[num] = val;
        }
        virtual std::vector<double> getParams()
        {
            return params;
        }
        //! Evaluate the force and energy
        virtual double getPairForceDivR(Eigen::Vector3d rij, double di, double dj)
        {
           return 0;
        }
        
        //! Evaluate the force and energy
        virtual double getPairEnergy(Eigen::Vector3d rij, double di, double dj)
        {
            return 0;
        }
        
        //! Evaluate the force and energy
        virtual double getBondStiffness(Eigen::Vector3d rij, double di, double dj)
        {
            return 0;
        }
        
        double scaled_rcut;
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

        double getRcut(Eigen::Vector3d rij, double di, double dj)
        {
            double rsq_ij = rij.dot(rij);
            double rsq_cut = scaled_rcut*scaled_rcut;
            ShortRangeModel model(rsq_ij, rsq_cut, params);
            return model.computeRcut(di,dj);
        };
        
        //! Evaluate the force and energy
        double getPairForceDivR(Eigen::Vector3d rij, double di, double dj)
        {
            double rsq_ij = rij.dot(rij);
            double rsq_cut = scaled_rcut*scaled_rcut;
            ShortRangeModel model(rsq_ij, rsq_cut, params);
            return model.computeForceDivR(di,dj);
        };
        
        //! Evaluate the force and energy
        double getPairEnergy(Eigen::Vector3d rij, double di, double dj)
        {
            double rsq_ij = rij.dot(rij);
            double rsq_cut = scaled_rcut*scaled_rcut;
            ShortRangeModel model(rsq_ij, rsq_cut, params);
            return model.computeEnergy(di,dj);
        };
        
        //! Evaluate the force and energy
        double getBondStiffness(Eigen::Vector3d rij, double di, double dj)
        {
            double rsq_ij = rij.dot(rij);
            double rsq_cut = scaled_rcut*scaled_rcut;
            ShortRangeModel model(rsq_ij, rsq_cut, params);
            return model.computeStiffness(di,dj);
        };
};

void export_PairPotential(py::module& m)
{
    py::class_<PairPotential, std::shared_ptr<PairPotential> > (m, "PairPotential")
    .def(py::init<std::vector<double> >())
    .def(py::init<double, std::vector<double> >())
    .def("getRcut", &PairPotential::getRcut)
    .def("getPairForceDivR", &PairPotential::getPairForceDivR)
    .def("getPairEnergy", &PairPotential::getPairEnergy)
    .def("getParams", &PairPotential::getParams)
    .def("setParams", &PairPotential::setParams)
    .def_readwrite("scaled_rcut", &PairPotential::scaled_rcut)
    ;
};

template < class T > 
void export_ShortRangePairPotential(py::module& m, const std::string& name)
{
    py::class_<T, PairPotential, std::shared_ptr<T> >(m, name.c_str())
    .def(py::init<std::vector<double> >())
    .def(py::init<double, std::vector<double> >())
    .def("getRcut", &T::getRcut)
    .def("getPairForceDivR", &T::getPairForceDivR)
    .def("getPairEnergy", &T::getPairEnergy)
    .def("getParams", &T::getParams)
    .def("setParams", &T::setParams)
    .def_readwrite("scaled_rcut", &T::scaled_rcut)
    ;
};

#endif // __PAIR_POTENTIAL_H__
