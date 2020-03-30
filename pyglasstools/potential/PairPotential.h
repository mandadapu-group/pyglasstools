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
            : m_di(0), m_dj(0), m_rij(Vector3d::Zero())
        {
        };
        
        //If parameter is given
        PairPotential(std::vector<double> _params)
            : m_di(0), m_dj(0), m_rij(Vector3d::Zero()), m_params(_params)
        {
        };
        virtual ~PairPotential(){};
        bool needsDiameter()
        { 
            return true; 
        };
        
        bool isPolydisperse()
        { 
            return true; 
        };

        void setDiameters(double di, double dj)
        {   
           m_di = di;
           m_dj = dj;
        };

        std::tuple<double, double> getDiameters()
        {
            return std::make_tuple(m_di,m_dj);  
        };

        void setRij(Vector3d rij)
        {
            m_rij = rij;
        };
        Vector3d getRij()
        {
            return m_rij;//m_rij = rij;
        };
        void setParams(std::vector<double> params)
        {
            m_params = params;
        };
        std::vector<double> getParams()
        {
            return m_params;
        };
        //Virtual methods here
        virtual void setScaledRcut(double rcut)
        {
        };
        virtual double getScaledRcut()
        {
            return 0.0;
        };
        virtual double getRcut()
        {
            return 0.0;
        };
        //! Evaluate the force and energy
        virtual double getPairForce()
        {
            return 0.0;
        };
    protected:
        double m_di;
        double m_dj;
        Vector3d m_rij;
        std::vector<double> m_params;
};

template<class ShortRangeModel>
class PYBIND11_EXPORT ShortRangePairPotential : public PairPotential
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        //Zero initialize everything
        ShortRangePairPotential(std::vector<double> _params)
            : PairPotential(_params), m_rcut(0)
        {
        };
        
        //Zero initialize everything
        ShortRangePairPotential(double _rcut, std::vector<double> _params)
            : PairPotential(_params), m_rcut(_rcut)
        {
        };

        ~ShortRangePairPotential(){};

        void setScaledRcut(double rcut)
        {
            m_rcut = rcut;
        };
        double getScaledRcut()
        {
            //py::print("scaled rcut is ", m_rcut);
            return m_rcut;
        };
        double getRcut()
        {
            double rsq_ij = m_rij.dot(m_rij);
            double rsq_cut = m_rcut*m_rcut;
            ShortRangeModel model(rsq_ij, rsq_cut, m_params);
            return model.computeRcut(m_di,m_dj);
        };
        //! Evaluate the force and energy
        double getPairForce()
        {
            double rsq_ij = m_rij.dot(m_rij);
            double rsq_cut = m_rcut*m_rcut;
            ShortRangeModel model(rsq_ij, rsq_cut, m_params);
            return model.computeForce(m_di,m_dj);
        };
    protected:
        double m_rcut;
};

void export_PairPotential(py::module& m)
{
    py::class_<PairPotential, std::shared_ptr<PairPotential> > (m, "PairPotential")
    .def(py::init<std::vector<double> >())
    .def("needsDiameter", &PairPotential::needsDiameter)
    .def("isPolydisperse", &PairPotential::isPolydisperse)
    .def("setDiameters", &PairPotential::setDiameters)
    .def("getDiameters", &PairPotential::getDiameters)
    .def("setRij", &PairPotential::setRij)
    .def("getRij", &PairPotential::getRij)
    .def("setParams", &PairPotential::setParams)
    .def("getParams", &PairPotential::getParams)
    .def("setScaledRcut", &PairPotential::setScaledRcut)
    .def("getScaledRcut", &PairPotential::getScaledRcut)
    .def("getRcut", &PairPotential::getRcut)
    .def("getPairForce", &PairPotential::getPairForce)
    ;
};

template < class T > 
void export_ShortRangePairPotential(py::module& m, const std::string& name)
{
    py::class_<T, PairPotential, std::shared_ptr<T> >(m, name.c_str())
    .def(py::init<double, std::vector<double> >())
    ;
};




#endif // __PAIR_POTENTIAL_H__
