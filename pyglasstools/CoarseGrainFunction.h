#ifndef __COARSE_GRAIN_FUNC_H__
#define __COARSE_GRAIN_FUNC_H__

#include "Quadrature.h"

#include <Eigen/Dense>
#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/eigen.h"
namespace py = pybind11;
using namespace Eigen;

//Encapsulates All Coarse-Graining Function Needs
template<class Distribution>
class PYBIND11_EXPORT CoarseGrainFunction
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        //Zero initialize 
        CoarseGrainFunction()  
            : m_rcut(0), m_x(Vector3d::Zero()), m_ri(Vector3d::Zero()), m_dr(Vector3d::Zero())
        {
        };
        //Parametrize initialization
        CoarseGrainFunction(double cg_rcut)  
            : m_rcut(cg_rcut), m_x(Vector3d::Zero()), m_ri(Vector3d::Zero()), m_dr(Vector3d::Zero())
        {
        };
        //Parametrize initialization
        CoarseGrainFunction(double cg_rcut, Vector3d x, Vector3d ri, Vector3d dr)  
            : m_rcut(cg_rcut), m_x(x), m_ri(ri), m_dr(dr)
        {
        };
        virtual ~CoarseGrainFunction(){};

        virtual void setRcut(double cg_rcut)
        {
            m_rcut = cg_rcut;
        };
        virtual double getRcut()
        {
            return m_rcut;
        };
        
        virtual void setX(Vector3d x)
        {
            m_x = x;
        };
        virtual Vector3d getX()
        {
            return m_x;
        };
        
        virtual void setRi(Vector3d ri)
        {
            m_ri = ri;
        };
        virtual Vector3d getRi()
        {
            return m_ri;
        };
        
        virtual void setDR(Vector3d dr)
        {
            m_dr = dr;
        };
        virtual Vector3d getDR()
        {
            return m_dr;
        };
        
        virtual double getDeltaFunc()
        {
            Vector3d dr = m_x-m_ri;
            double dr_sq = dr.dot(dr);
            Distribution func(dr_sq, m_rcut);
            return func.compute();
        };
        
        virtual double getObjFunc(double s)
        {
            Vector3d dr = m_x-(m_ri+s*m_dr);
            double dr_sq = dr.dot(dr);
            Distribution func(dr_sq, m_rcut);
            return func.compute();
        };

        virtual double getBondFunc()
        {
            return GSLQuadrature([&](double s) { return getObjFunc(s); }, {0,1});
        };

    private:
        double m_rcut;
        Vector3d m_x;
        Vector3d m_ri;
        Vector3d m_dr;
};

template < class T > void export_CoarseGrainFunction(py::module& m, const std::string& name)
{
    py::class_<T, std::shared_ptr<T> >(m, name.c_str())
    .def(py::init<double>())
    .def(py::init<double, Vector3d, Vector3d, Vector3d>())
    .def("setRcut", &T::setRcut)
    .def("getRcut", &T::getRcut)
    .def("setX", &T::setX)
    .def("getX", &T::getX)
    .def("setRi", &T::setRi)
    .def("getRi", &T::getRi)
    .def("setDR", &T::setDR)
    .def("getDR", &T::getDR)
    .def("getDeltaFunc", &T::getDeltaFunc)
    .def("getObjFunc", &T::getObjFunc)
    .def("getBondFunc", &T::getBondFunc)
    ;
};

//(1): Coarse-Graining Function based on an 8-th order polynomial
class Octic
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        Octic(double _dr_sq, double _rcut)
            : dr_sq(_dr_sq), rcut(_rcut)
            {
            }
        ~Octic(){}; 
        //! Evaluate the force and energy
        virtual double compute()
        {
            double rcutsq = rcut*rcut;
            
            // compute the force divided by r in force_divr
            if (dr_sq < rcutsq)
            {
                double r4 = dr_sq*dr_sq/(rcutsq*rcutsq);
                double r8 = r4*r4;
                return 15.0/(8*M_PI*rcutsq)*(1-2*r4+r8);
            }
            else
                return 0.0;
        }

    protected:
        double dr_sq;     //!< Stored rsq from the constructor
        double rcut;  //!< Stored rcutsq from the constructor
};
#endif //__COARSE_GRAIN_FUNC_H__
