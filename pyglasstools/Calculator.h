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
#include <pyglasstools/observables/Observables.h>
#include <pyglasstools/cgfunc/CoarseGrainFunction.h>

#include "extern/pybind11/include/pybind11/pybind11.h"
namespace py = pybind11;

#include <Aboria.h>
namespace abr = Aboria;

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Core>

class PYBIND11_EXPORT GlobalCalculator
{
    public:
        GlobalCalculator(std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential)
            : m_sysdata(sysdata), m_potential(potential)
            {
                double maxforce_rcut = potential->scaled_rcut;
                double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                                        std::end(abr::get<diameter>(m_sysdata->particles)) );
                maxforce_rcut *= maxdiameter;
                max_rcut = std::max(maxforce_rcut,maxdiameter); //whichever is maximum
            };
        ~GlobalCalculator(){};

        virtual void addObservable(const std::shared_ptr<Observable>& obs)
            {
                m_observables.insert(std::pair<std::string, std::shared_ptr<Observable> >(obs->name, obs));
            }
        virtual MatrixXd getGlobalObservable(const std::string& name)
            {
                return m_observables[name]->getGlobalValue();
            }
            
        virtual void compute();
        
        virtual void computeLocalObs(const AboriaParticles::value_type& particle_i) 
            {
                for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
                {
                    if (it->second->islocal && !it->second->isfield)
                        it->second->accumulate(particle_i);
                    else 
                        continue;
                }
            }
        virtual void computePairObs(const AboriaParticles::value_type& particle_i, 
                                    const AboriaParticles::value_type& particle_j)
            {
                for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
                {
                    if (it->second->useforce && !it->second->islocal && !it->second->isfield)
                        it->second->accumulate(particle_i,particle_j,m_potential);
                    else if (!it->second->useforce && !it->second->islocal && !it->second->isfield)
                        it->second->accumulate(particle_i,particle_j);
                }
            }
        virtual void clearState()
            {
                for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
                {
                    it->second->clear();
                }
            }
    protected:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< PairPotential > m_potential;
        double max_rcut;
        std::map< std::string, std::shared_ptr<Observable> > m_observables;

};

//Compute a Global Observable
void GlobalCalculator::compute()
{
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
                m_potential->rij = rij;
                
                //Don't forget to set diameters of the potential
                m_potential->di =  abr::get<diameter>(*p_i);
                m_potential->dj =  abr::get<diameter>(*p_j);
                
                computePairObs(*p_i,*p_j);
            }
        }
    } //end of for loop n
};

class PYBIND11_EXPORT LocalCalculator : public GlobalCalculator
{
    public:
        LocalCalculator(    std::shared_ptr< ParticleSystem > sysdata,
                            std::shared_ptr< PairPotential > potential, 
                            std::shared_ptr< CoarseGrainFunction > cgfunc )
            : GlobalCalculator(sysdata,potential), m_cgfunc(cgfunc)
        {
        }; 
        ~LocalCalculator(){};
        
        virtual void computelocal(const std::vector< Eigen::Vector3d >& gridpoints);
        
        std::vector< MatrixXd > getField(std::string name)
            {
                return m_observables[name]->getField();
            }

        virtual void computeLocalObsPerGrid(AboriaParticles::value_type particle_i, double cgval, unsigned int grid_id) 
            {
                for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
                {
                    if (it->second->islocal && it->second->isfield)
                        it->second->accumulate(particle_i,cgval,grid_id);
                    else 
                        continue;
                }
            }
        virtual void computePairObsPerGrid( AboriaParticles::value_type particle_i, 
                                            AboriaParticles::value_type particle_j, double bondval, unsigned int grid_id)
            {
                for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
                {
                    if (it->second->useforce && !it->second->islocal && it->second->isfield)
                        it->second->accumulate(particle_i,particle_j,m_potential, bondval, grid_id);
                    else if (!it->second->useforce && !it->second->islocal && it->second->isfield)
                        it->second->accumulate(particle_i,particle_j,bondval, grid_id);
                }
            }
    private:
        std::shared_ptr< CoarseGrainFunction > m_cgfunc; //!< particle system, equipped with neighbor list 
};

//Compute a Global Observable
void LocalCalculator::computelocal(const std::vector< Eigen::Vector3d >& gridpoints)
{
    clearState();
    for (unsigned int i = 0; i < gridpoints.size(); ++i)
    {
        for( auto p_i = abr::euclidean_search(m_sysdata->particles.get_query(), 
                        abr::vdouble3(gridpoints[i][0],gridpoints[i][1],gridpoints[i][2]), m_cgfunc->cg_rcut); 
                        p_i != false; ++p_i)
        {
            
            //Set grid point X and position of particle Ri
            m_cgfunc->x = gridpoints[i];
            Vector3d dr(p_i.dx()[0], p_i.dx()[1],p_i.dx()[2]);
            m_cgfunc->ri = gridpoints[i]-dr;
            double cgval = m_cgfunc->getDeltaFunc();
            computeLocalObsPerGrid(*p_i,cgval, i);
            
            //Compute a list of local obsercavles 
            //Next we loop through the j-th particles for the virial stress
            for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                            abr::get<position>(*p_i), max_rcut); 
                            p_j != false; ++p_j)
            {
                //Make sure the particle is unique
                if (abr::get<abr::id>(*p_i) != abr::get<abr::id>(*p_j))
                {
                    //set the distance between particle i and particle j
                    Eigen::Vector3d _rij(-p_j.dx()[0], -p_j.dx()[1], -p_j.dx()[2]);
                    m_potential->rij = _rij;
                    m_cgfunc->rij = _rij;
                    double bondval = m_cgfunc->getBondFunc();
                    
                    //Don't forget to set diameters of the potential
                    m_potential->di =  abr::get<diameter>(*p_i);
                    m_potential->dj =  abr::get<diameter>(*p_j);
                    computePairObsPerGrid(*p_i,*p_j,bondval, i);
                } 
            } 
        }
    }
};

void export_GlobalCalculator(py::module& m)
{
    py::class_<GlobalCalculator, std::shared_ptr<GlobalCalculator> >(m,"GlobalCalculator")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >  >())
    .def("addObservable", &GlobalCalculator::addObservable)
    .def("getGlobalObservable", &GlobalCalculator::getGlobalObservable)
    .def("compute", &GlobalCalculator::compute)
    ;
};

void export_LocalCalculator(py::module& m)
{
    py::class_<LocalCalculator, GlobalCalculator, std::shared_ptr<LocalCalculator> >(m,"LocalCalculator")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< CoarseGrainFunction > >())
    .def("getField", &LocalCalculator::getField)
    .def("computelocal", &LocalCalculator::computelocal)
    ;
};

#endif //__CALCULATOR_H__
