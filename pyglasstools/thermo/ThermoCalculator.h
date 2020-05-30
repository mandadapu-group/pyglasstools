#ifndef __THERMO_CALCULATOR_H__
#define __THERMO_CALCULATOR_H__

#include <pyglasstools/Calculator.h>
#include "ThermoProperty.h"

class PYBIND11_EXPORT ThermoCalculator : public Calculator
{
    private:
        std::map< std::string, std::shared_ptr< ThermoProperty > > m_observables;
    public:

        ThermoCalculator(   std::shared_ptr< ParticleSystem > sysdata, 
                            std::shared_ptr< PairPotential > potential)
            : Calculator(sysdata, potential){};
        ~ThermoCalculator(){};

        virtual void addObservable(std::shared_ptr< ThermoProperty > obs)
        {
            m_observables.insert(std::pair<std::string, std::shared_ptr< ThermoProperty > >(obs->name, obs));
        }
        virtual void compute();
        
        virtual void computeLocalObs(const AboriaParticles::value_type& particle_i) 
        {
            for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
            {
                if (it->second->islocal)
                    it->second->accumulate(particle_i);
                else 
                    continue;
            }
        }
        
        virtual void computePairObs(const AboriaParticles::value_type& particle_i, 
                                    const AboriaParticles::value_type& particle_j,
                                    Eigen::Vector3d rij)
        {
            for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
            {
                if (!it->second->islocal)
                {
                    it->second->accumulate(particle_i,particle_j,rij,m_potential);
                }
                else
                    continue;
            }
        }
        
        virtual void clearState()
        {
            for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
            {
                it->second->clear();
            }
        }
        
        void divideByVolume()
        {
            for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
            {
                it->second->divideByVolume(m_sysdata->simbox->vol);
            }
        }
};

//Compute a Global Observable
void ThermoCalculator::compute()
{
    clearState();
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
                Eigen::Vector3d rij(p_j.dx()[0], p_j.dx()[1], p_j.dx()[2]);
                computePairObs(*p_i,*p_j,rij);
            }
        } // end of for loop on j
    } //end of for loop on i
    divideByVolume();
};

void export_ThermoCalculator(py::module& m)
{
    py::class_<ThermoCalculator, Calculator, std::shared_ptr<ThermoCalculator> >(m,"ThermoCalculator")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >  >())
    .def("compute", &ThermoCalculator::compute)
    .def("setSystemData", &ThermoCalculator::setSystemData)
    .def("addObservable", &ThermoCalculator::addObservable)
    ;
};

#endif //__THERMO_CALCULATOR_H__
