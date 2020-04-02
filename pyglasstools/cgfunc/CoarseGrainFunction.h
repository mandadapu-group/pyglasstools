#ifndef __COARSE_GRAIN_FUNC_H__
#define __COARSE_GRAIN_FUNC_H__

#include "Quadrature.h"

#include "../extern/pybind11/include/pybind11/pybind11.h"
#include "../extern/pybind11/include/pybind11/eigen.h"
namespace py = pybind11;

#include <Eigen/Dense>

//Encapsulates All Coarse-Graining Function Needs
class PYBIND11_EXPORT CoarseGrainFunction
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        //Zero initialize 
        CoarseGrainFunction()  
            : cg_rcut(0), x(Eigen::Vector3d::Zero()), ri(Eigen::Vector3d::Zero()), rij(Eigen::Vector3d::Zero())
        {
        };
        CoarseGrainFunction(double _cg_rcut)  
            : cg_rcut(_cg_rcut), x(Eigen::Vector3d::Zero()), ri(Eigen::Vector3d::Zero()), rij(Eigen::Vector3d::Zero())
        {
        };
        //Parametrize initialization
        CoarseGrainFunction(double _cg_rcut, Eigen::Vector3d _x, Eigen::Vector3d _ri, Eigen::Vector3d _rij)  
            : cg_rcut(_cg_rcut), x(_x), ri(_ri), rij(_rij)
        {
        };
        virtual ~CoarseGrainFunction(){};
        
        virtual double getDeltaFunc()
        {
            return 0.0;
        };
        
        virtual double getObjFunc(double s)
        {
            return 0.0;
        };

        virtual double getBondFunc()
        {
            return 0.0;
        };
        
        double cg_rcut;
        Eigen::Vector3d x;
        Eigen::Vector3d ri;
        Eigen::Vector3d rij;
};

template<class Distribution>
class PYBIND11_EXPORT ShortRangeCGFunc : public CoarseGrainFunction
{
    public:
        //Zero initialize 
        ShortRangeCGFunc(){};

        //Parametrize initialization
        ShortRangeCGFunc(double cg_rcut)  
            : CoarseGrainFunction(cg_rcut)
        {
        };
        //Parametrize initialization
        ShortRangeCGFunc(double cg_rcut, Eigen::Vector3d x, Eigen::Vector3d ri, Eigen::Vector3d rij)  
            : CoarseGrainFunction(cg_rcut, x,ri,rij)
        {
        };
        double getDeltaFunc()
        {
            Eigen::Vector3d dr = x-ri;
            double dr_sq = dr.dot(dr);
            Distribution func(dr_sq, cg_rcut);
            return func.compute();
        };
        
        double getObjFunc(double s)
        {
            Eigen::Vector3d dr = x-(ri+s*rij);
            double dr_sq = dr.dot(dr);
            Distribution func(dr_sq, cg_rcut);
            return func.compute();
        };

        double getBondFunc()
        {
            return GSLQuadrature([&](double s) { return getObjFunc(s); }, {0,1});
        };
};

void export_CoarseGrainFunction(py::module& m)
{
    py::class_<CoarseGrainFunction, std::shared_ptr<CoarseGrainFunction> >(m, "CoarseGrainFunction")
    .def(py::init<double>())
    .def(py::init<double, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>())
    .def("getDeltaFunc", &CoarseGrainFunction::getDeltaFunc)
    .def("getObjFunc", &CoarseGrainFunction::getObjFunc)
    .def("getBondFunc", &CoarseGrainFunction::getBondFunc)
    .def_readwrite("cg_rcut", &CoarseGrainFunction::cg_rcut)
    .def_readwrite("x", &CoarseGrainFunction::x)
    .def_readwrite("ri", &CoarseGrainFunction::ri)
    .def_readwrite("rij", &CoarseGrainFunction::rij)
    ;
};

template < class T > 
void export_ShortRangeCGFunc(py::module& m, const std::string& name)
{
    py::class_<T, CoarseGrainFunction, std::shared_ptr<T> >(m, name.c_str())
    .def(py::init<double>())
    .def(py::init<double, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>())
    ;
};
#endif //__COARSE_GRAIN_FUNC_H__
