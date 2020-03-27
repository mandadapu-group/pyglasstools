#ifndef __PAIR_POTENTIAL_H__
#define __PAIR_POTENTIAL_H__

#include <Eigen/Dense>
#include <cmath>
#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/eigen.h"
#include "extern/pybind11/include/pybind11/stl.h"
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

/*
//Basic interface for common functions
class PYBIND11_EXPORT PyPairPotential : public PairPotential
{
    public:
        //Zero initialize everything
        using PairPotential::PairPotential;
        
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
};
*/

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
    //.def("setScaledRcut", &T::setScaledRcut)
    //.def("getScaledRcut", &T::getScaledRcut)
   // .def("getRcut", &T::getRcut)
   // .def("getPairForce", &T::getPairForce)
    ;
};

// DEFINE Your Models HERE 
//(1): 1-parameter LennardJones
class LennardJones
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        LennardJones(double _rsq, double _rcutsq, std::vector<double> _params)
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
        virtual double computeRcut(double d_i, double d_j)
        {
            double sigma = 0.5*(d_i+d_j);
            return rcutsq*sigma*sigma;
        }

    protected:
        double rsq;     //!< Stored rsq from the constructor
        double rcutsq;  //!< Stored rcutsq from the constructor
        double eps;     //!< epsilon_parameter
};

class Polydisperse12
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        Polydisperse12(double _rsq, double _rcutsq, std::vector<double> _params)
            : rsq(_rsq), rcutsq(_rcutsq), v0(_params[0]), eps(_params[1])
            {
                c0 =  -28.0*v0/pow(_rcutsq,6);
                c1 =  48.0*v0/pow(_rcutsq,7);
                c2 =  -21.0*v0/pow(_rcutsq,8);
            }
        ~Polydisperse12(){}; 
        //! Evaluate the force and energy
        virtual double computeForce(double d_i, double d_j)
        {
            double sigma = 0.5*(d_i+d_j)*(1-eps*fabs(d_i-d_j));
            double actualrcutsq = rcutsq*sigma*sigma;
            
            // compute the force divided by r in force_divr
            if (rsq < actualrcutsq)
                {
                    double r2inv = 1.0/rsq;
                    r2inv *= sigma*sigma;
                    double _rsq = rsq/(sigma*sigma);
                    double r6inv = r2inv * r2inv * r2inv;
                    double force_divr = 12.0*v0*r2inv*r6inv*r6inv-2.0*c1 -4.0*c2*_rsq;
                    return force_divr;
                }
            else
                return 0.0;
        }
        virtual double computeRcut(double d_i, double d_j)
        {
            double sigma = 0.5*(d_i+d_j)*(1-eps*fabs(d_i-d_j));
            return rcutsq*sigma*sigma;
        }

    protected:
        double rsq;     //!< Stored rsq from the constructor
        double rcutsq;  //!< Stored rcutsq from the constructor
        double v0;
        double eps;     //!< epsilon_parameter
        double c0;
        double c1;
        double c2;
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
        ForceShiftedLennardJones(double _rsq, double _rcutsq, std::vector<double> _params)
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
