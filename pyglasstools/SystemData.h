#ifndef __SYSTEM_DATA_H__
#define __SYSTEM_DATA_H__

//#include "extern/aabbcc/src/AABB.h"
#include <Eigen/Dense>
#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/eigen.h"
namespace py = pybind11;
using namespace Eigen;

class PYBIND11_EXPORT SystemData
{
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        SystemData(VectorXd diameter, VectorXd mass, Matrix<double, Dynamic, 3, RowMajor> atomposition, Matrix<double, Dynamic, 3, RowMajor> atomvelocity)
        {
            m_mass = mass;
            m_diameter = diameter;  
            m_atomposition = atomposition;
            m_atomvelocity = atomvelocity;
                  
        };
        ~SystemData(){};
        void setMass(VectorXd mass)
        {
            int rowsize = mass.rows();
            m_mass.resize(rowsize);
            m_mass = mass; 
        };
        VectorXd getMass()
        {
            return m_mass; 
        };
        
        void setDiameter(VectorXd diameter)
        {
            int rowsize = diameter.rows();
            m_diameter.resize(rowsize);
            m_diameter = diameter; 
        };
        VectorXd getDiameter()
        {
            return m_diameter; 
        };
        
        void setAtomPosition(Matrix<double, Dynamic, 3, RowMajor> atomposition)
        {
            int rowsize = atomposition.rows();
            m_atomposition.resize(rowsize,3);
            m_atomposition = atomposition; 
        };
        Matrix<double, Dynamic, 3, RowMajor> getAtomPosition()
        {
            return m_atomposition; 
        };
        
        void setAtomVelocity(Matrix<double, Dynamic, 3, RowMajor> atomvelocity)
        {
            int rowsize = atomvelocity.rows();
            m_atomvelocity.resize(rowsize,3);
            m_atomvelocity = atomvelocity; 
        };
        Matrix<double, Dynamic, 3, RowMajor> getAtomVelocity()
        {
            return m_atomvelocity; 
        };

    private: 
        //Atomic properties, fed in to the system
        VectorXd m_mass;
        VectorXd m_diameter;
        Matrix<double, Dynamic, 3, RowMajor> m_atomposition;
        Matrix<double, Dynamic, 3, RowMajor> m_atomvelocity;
};

// I might need to do something special since I'm using Eigen types
//an export function here
void export_SystemData(py::module& m)
{
    py::class_<SystemData>(m,"SystemData")
    .def(py::init<VectorXd, VectorXd, Matrix<double, Dynamic, 3, RowMajor>, Matrix<double, Dynamic, 3, RowMajor>>())
    .def("getMass", &SystemData::getMass)
    .def("setMass", &SystemData::setMass)
    .def("getDiameter", &SystemData::getDiameter)
    .def("setDiameter", &SystemData::setDiameter)
    .def("getAtomPosition", &SystemData::getAtomPosition)
    .def("setAtomPosition", &SystemData::setAtomPosition)
    .def("getAtomVelocity", &SystemData::getAtomVelocity)
    .def("setAtomVelocity", &SystemData::setAtomVelocity)
    ;
};
#endif //__SYSTEM_DATA_H__
