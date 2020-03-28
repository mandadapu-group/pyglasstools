#ifndef __CALCULATOR_H__
#define __CALCULATOR_H__
#include <vector>
#include <cmath>
#include <string>
#include <map>

#include <algorithm>
#include <memory>

#include "MathAndTypes.h"
#include "ParticleSystem.h"
#include <pyglasstools/cgfunc/CoarseGrainFunction.h>
#include "SimBox.h"
#include <pyglasstools/observables/Observables.h>

#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/stl.h"
#include "extern/pybind11/include/pybind11/eigen.h"
namespace py = pybind11;

#include <Aboria.h>
namespace abr = Aboria;

#include <Eigen/Dense>
#include <Eigen/StdVector>

/*
class PYBIND11_EXPORT GridPoints
{
    public:
        GridPoints() : totsize(0){};
        GridPoints(const std::vector< Eigen::Vector3d >& _points)
            : totsize(_points.size()), points(_points)
        {};
        ~GridPoints(){};
        
        unsigned int totsize;
        std::vector< Eigen::Vector3d > points;
};

void export_GridPoints(py::module& m)
{
    py::class_<GridPoints, std::shared_ptr<GridPoints> >(m,"GridPoints")
    .def(py::init<std::vector< Eigen::Vector3d > >())
    ;
};
*/

class PYBIND11_EXPORT Calculator
{
    public:
        Calculator(std::shared_ptr< ParticleSystem > sysdata) : m_sysdata(sysdata){};
        ~Calculator(){};
        
        void addObservable(std::shared_ptr<Observable> obs)
            {
                m_observables.insert(std::pair<std::string, std::shared_ptr<Observable> >(obs->name, obs));
            }
        
        MatrixXd getObservable(std::string name)
            {
                return m_observables[name]->val;
            }
            
        virtual void compute(){};
        
    protected:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::map< std::string, std::shared_ptr<Observable> > m_observables; 
};


class PYBIND11_EXPORT GlobalCalculator : public Calculator
{
    public:
        GlobalCalculator(std::shared_ptr< ParticleSystem > sysdata)
            : Calculator(sysdata)
        {
            double maxforce_rcut = m_sysdata->potential->getScaledRcut();
            double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                                    std::end(abr::get<diameter>(m_sysdata->particles)) );
            maxforce_rcut *= maxdiameter;
            max_rcut = std::max(maxforce_rcut,maxdiameter); //whichever is maximum
        };
        ~GlobalCalculator(){};
        void compute();
        void computeLocalObs(AboriaParticles::value_type particle_i) 
            {
                for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
                {
                    if (it->second->islocal)
                        it->second->accumulate(particle_i);
                    else 
                        continue;
                }
            }
        void computePairObs(    AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j ) 
            {
                for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
                {
                    if (it->second->islocal)
                        continue;
                    else
                    {
                        if (it->second->useforce)
                        {
                            it->second->accumulate(particle_i,particle_j,m_sysdata->potential);
                        }
                        else
                            it->second->accumulate(particle_i,particle_j);
                    }
                }
            }
    private:
        double max_rcut;
};

//Compute a Global Observable
void GlobalCalculator::compute()
{
    //#pragma omp parallel for reduction(+:totTvxy)
    for( auto p_i = m_sysdata->particles.begin(); p_i != m_sysdata->particles.end(); ++p_i)
    {
        //Compute a list of local obsercavles 
        computeLocalObs(*p_i);

        //Next we loop through the j-th particles for the virial stress
        for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                    abr::get<position>(*p_i), max_rcut); p_j != false; ++p_j)
        {
            //Make sure the particle is unique
            if (abr::get<abr::id>(*p_i) < abr::get<abr::id>(*p_j))
            {
                //set the distance between particle i and particle j
                Eigen::Vector3d rij;
                rij << p_j.dx()[0], p_j.dx()[1], p_j.dx()[2];
                m_sysdata->potential->setRij(rij);
                
                //Don't forget to set diameters of the potential
                m_sysdata->potential->setDiameters( abr::get<diameter>(*p_i),
                                         abr::get<diameter>(*p_j));
                computePairObs(*p_i,*p_j);
            }
        }
    } //end of for loop n
};

void export_Calculator(py::module& m)
{
    py::class_<Calculator, std::shared_ptr<Calculator> >(m,"Calculator")
    .def(py::init< std::shared_ptr< ParticleSystem > >())
    .def("addObservable", &Calculator::addObservable)
    .def("getObservable", &Calculator::getObservable)
    .def("compute", &Calculator::compute)
    ;
};

void export_GlobalCalculator(py::module& m)
{
    py::class_<GlobalCalculator, Calculator, std::shared_ptr<GlobalCalculator> >(m,"GlobalCalculator")
    .def(py::init< std::shared_ptr< ParticleSystem > >())
    ;
};

#endif //__CALCULATOR_H__
