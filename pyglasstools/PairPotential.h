#ifndef __PAIR_POTENTIAL_H__
#define __PAIR_POTENTIAL_H__

#include <Eigen/Dense>
#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/eigen.h"
namespace py = pybind11;
using namespace Eigen;

//Basic interface for common functions
//Forward declare parameters struct
template<class Model>
class PYBIND11_EXPORT PairPotential
{
    public:
        //Zero initialize everything
        PairPotential(VectorXd _params)
            : m_di(0), m_dj(0), m_rij(Vector3d::Zero()), m_rcut(0), m_params(_params)
        {
        };
        
        //Zero initialize everything
        PairPotential(double _rcut, VectorXd _params)
            : m_di(0), m_dj(0), m_rij(Vector3d::Zero()), m_rcut(_rcut), m_params(_params)
        {
        };

        ~PairPotential(){};
        
        virtual bool needsDiameter()
        { 
            return true; 
        };
        
        virtual bool isPolydisperse()
        { 
            return true; 
        };

        virtual void setDiameters(double di, double dj)
        {   
           m_di = di;
           m_dj = dj;
        };

        virtual std::tuple<double, double> getDiameters()
        {
            return std::make_tuple(m_di,m_dj);  
        };

        virtual void setRij(Vector3d rij)
        {
            m_rij = rij;
        };
        virtual Vector3d getRij()
        {
            return m_rij;//m_rij = rij;
        };

        virtual void setRcut(double rcut)
        {
            m_rcut = rcut;
        };
        virtual double getRcut()
        {
            return m_rcut;
        };
        
        virtual void setParams(VectorXd params)
        {
            m_params = params;
        };
        virtual VectorXd getParams()
        {
            return m_params;
        };

        //! Evaluate the force and energy
        virtual double getPairForce()
        {
            double rsq_ij = m_rij.dot(m_rij);
            double rsq_cut = m_rcut*m_rcut;
            Model model(rsq_ij, rsq_cut, m_params);
            return model.computeForce(m_di,m_dj);
        };
    protected:
        double m_di;
        double m_dj;
        Vector3d m_rij;
        double m_rcut;
        VectorXd m_params;
};

template < class T > void export_PairPotential(py::module& m, const std::string& name)
{
    py::class_<T, std::shared_ptr<T> >(m, name.c_str())
    .def(py::init<VectorXd>())
    .def(py::init<double, VectorXd>())
    .def("needsDiameter", &T::needsDiameter)
    .def("isPolydisperse", &T::isPolydisperse)
    .def("setDiameters", &T::setDiameters)
    .def("getDiameters", &T::getDiameters)
    .def("setRij", &T::setRij)
    .def("getRij", &T::getRij)
    .def("setRcut", &T::setRcut)
    .def("getRcut", &T::getRcut)
    .def("setParams", &T::setParams)
    .def("getParams", &T::getParams)
    .def("getPairForce", &T::getPairForce)
    ;
};

//(1): 1-parameter LennardJones
class LennardJones
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        LennardJones(double _rsq, double _rcutsq, VectorXd _params)
            : rsq(_rsq), rcutsq(_rcutsq), eps(_params[0])
            {
            }
        ~LennardJones(){}; 
        //! Evaluate the force and energy
        virtual double computeForce(double d_i, double d_j)
        {
            double sigma = 0.5*(d_i+d_j);
            double actualrcutsq = rcutsq*sigma*sigma;
            
            // compute the force divided by r in force_divr
            if (rsq < actualrcutsq && eps != 0)
                {
                double r2inv = 1.0/rsq;
                r2inv *= sigma*sigma;
                double r6inv = r2inv * r2inv * r2inv;
                double force_divr = 4*eps*r2inv * r6inv * (12.0*r6inv - 6.0);
                return force_divr;
                }
            else
                return 0.0;
        }

    protected:
        double rsq;     //!< Stored rsq from the constructor
        double rcutsq;  //!< Stored rcutsq from the constructor
        double eps;     //!< epsilon_parameter
};

//(2): 1-parameter Force-Shifted LennardJones
class ForceShiftedLennardJones : public LennardJones
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        ForceShiftedLennardJones(double _rsq, double _rcutsq, VectorXd _params)
            : LennardJones(_rsq, _rcutsq, _params)
            {
            }
        ~ForceShiftedLennardJones(){};        
        //! Evaluate the force and energy
        virtual double computeForce(double d_i, double d_j)
        {
            double sigma = 0.5*(d_i+d_j);
            double actualrcutsq = rcutsq*sigma*sigma;
            
            // compute the force divided by r in force_divr
            if (rsq < actualrcutsq && eps != 0)
            {
                double r2inv = 1.0/rsq;
                r2inv *= sigma*sigma;
                double r6inv = r2inv * r2inv * r2inv;
                double force_divr = 4*eps*r2inv * r6inv * (12.0*r6inv - 6.0);
                r2inv = 1.0/rcutsq;
                r6inv = r2inv * r2inv * r2inv;
                force_divr -= 4*eps*r2inv * r6inv * (12.0*r6inv - 6.0);
                return force_divr;
            }
            else
                return 0.0;
        }
};

#endif // __PAIR_POTENTIAL_H__
