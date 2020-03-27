#ifndef __PAIR_POTENTIAL_H__
#define __PAIR_POTENTIAL_H__

#include <Eigen/Dense>
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
        };
        virtual double getRcut()
        {
            return 0.0;
        };
        virtual bool checkRange(Vector3d dr)
        {
            return true;
        };

        virtual void setParams(std::vector<double> params)
        {
            m_params = params;
        };
        virtual std::vector<double> getParams()
        {
            return m_params;
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

        virtual void setRcut(double rcut)
        {
            m_rcut = rcut;
        };
        virtual double getRcut()
        {
            return m_rcut;
        };
        virtual bool checkRange(Vector3d dr)
        {
            if (dr.dot(dr) < 0.25*m_rcut*m_rcut*(m_di+m_dj))
                return true;
            else
                return false;
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
    .def("getPairForce", &PairPotential::getPairForce)
    ;
};

template < class T > 
void export_ShortRangePairPotential(py::module& m, const std::string& name)
{
    py::class_<T, PairPotential, std::shared_ptr<T> >(m, name.c_str())
    .def(py::init<double, std::vector<double> >())
    .def("setRcut", &T::setRcut)
    .def("getRcut", &T::getRcut)
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
