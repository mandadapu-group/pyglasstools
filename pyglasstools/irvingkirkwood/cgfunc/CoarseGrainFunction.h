#ifndef __COARSE_GRAIN_FUNC_H__
#define __COARSE_GRAIN_FUNC_H__

#include "FixedPointQuadrature.h"
#include "Quadrature.h"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <Eigen/Dense>

//Encapsulates All Coarse-Graining Function Needs
class PYBIND11_EXPORT CoarseGrainFunction
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        //Zero initialize 
        CoarseGrainFunction()  
            : cg_rcut(0)
        {
        };
        CoarseGrainFunction(double _cg_rcut)  
            : cg_rcut(_cg_rcut)
        {
        };
        virtual ~CoarseGrainFunction(){};
        
        virtual void setRcut(double rcut)
        {
            cg_rcut = rcut;
        } 
        virtual double getRcut()
        {
            return cg_rcut;
        }
        virtual double getDeltaFunc(Eigen::Vector3d x, Eigen::Vector3d ri)
        {
            return 0.0;
        };
        
        double getObjFunc(double s, Eigen::Vector3d x, Eigen::Vector3d ri, Eigen::Vector3d rij)
        {
            return getDeltaFunc(x,ri+s*rij);
        };

        virtual double getBondFunc(Eigen::Vector3d x, Eigen::Vector3d ri, Eigen::Vector3d rij)
        {
            return 0.0;
        };
    
    protected:
        double cg_rcut;
};

template<class Distribution>
class PYBIND11_EXPORT FixedPointCGFunc : public CoarseGrainFunction
{
    public:
        //Zero initialize 
        FixedPointCGFunc(): quadorder(1) {};

        //Parametrize initialization
        FixedPointCGFunc(int order, double cg_rcut)  
            : CoarseGrainFunction(cg_rcut), quadorder(order), func(cg_rcut)
        {
        };

        virtual void setRcut(double rcut)
        {
            cg_rcut = rcut;
            func.setRcut(cg_rcut);
        } 
        virtual double getRcut()
        {
            return cg_rcut;
        }
        double getDeltaFunc(Eigen::Vector3d x, Eigen::Vector3d ri)
        {
            Eigen::Vector3d dr = x-ri;
            double dr_sq = dr.dot(dr);
            return func.compute(dr_sq);
        };
        
        double getBondFunc(Eigen::Vector3d x, Eigen::Vector3d ri, Eigen::Vector3d rij)
        {
            return GSLFixedPointQuadrature([&](double s) { return getObjFunc(s,x,ri,rij); }, {0,1}, quadorder);
        };
    private:
        unsigned int quadorder;
        Distribution func;
};

template<class Distribution>
class PYBIND11_EXPORT AdaptiveCGFunc : public CoarseGrainFunction
{
    public:
        //Zero initialize 
        AdaptiveCGFunc(): neval(1), epsrelerr(1e-9), epsabserr(1e-9) {};

        //Parametrize initialization
        AdaptiveCGFunc(int neval, double epsrelerr, double epsabserr, double cg_rcut)  
            : CoarseGrainFunction(cg_rcut), neval(neval), epsrelerr(epsrelerr), epsabserr(epsabserr), func(cg_rcut)
        {
        };
        virtual void setRcut(double rcut)
        {
            cg_rcut = rcut;
            func.setRcut(cg_rcut);
        } 
        virtual double getRcut()
        {
            return cg_rcut;
        }
        
        double getDeltaFunc(Eigen::Vector3d x, Eigen::Vector3d ri)
        {
            Eigen::Vector3d dr = x-ri;
            double dr_sq = dr.dot(dr);
            return func.compute(dr_sq);
        };
        

        double getBondFunc(Eigen::Vector3d x, Eigen::Vector3d ri, Eigen::Vector3d rij)
        {
            return GSLQuadrature([&](double s) { return getObjFunc(s,x,ri,rij); }, {0,1}, epsabserr, epsrelerr, neval);
        };
        
    private:
        unsigned int neval;
        double epsrelerr; 
        double epsabserr; 
        Distribution func;
};

void export_CoarseGrainFunction(py::module& m)
{
    py::class_<CoarseGrainFunction, std::shared_ptr<CoarseGrainFunction> >(m, "CoarseGrainFunction")
    .def(py::init<double>())
    ;
};

template < class T > 
void export_FixedPointCGFunc(py::module& m, const std::string& name)
{
    py::class_<T, CoarseGrainFunction, std::shared_ptr<T> >(m, name.c_str())
    .def(py::init<int, double>())
    .def("getRcut", &T::getRcut)
    .def("setRcut", &T::setRcut)
    .def("getDeltaFunc", &T::getDeltaFunc)
    .def("getObjFunc", &T::getObjFunc)
    .def("getBondFunc", &T::getBondFunc)
    ;
};

template < class T > 
void export_AdaptiveCGFunc(py::module& m, const std::string& name)
{
    py::class_<T, CoarseGrainFunction, std::shared_ptr<T> >(m, name.c_str())
    .def(py::init<int, double, double, double>())
    .def("getRcut", &T::getRcut)
    .def("setRcut", &T::setRcut)
    .def("getDeltaFunc", &T::getDeltaFunc)
    .def("getObjFunc", &T::getObjFunc)
    .def("getBondFunc", &T::getBondFunc)
    ;
};

#endif //__COARSE_GRAIN_FUNC_H__
