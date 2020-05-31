#ifndef __IRVING_KIRKWOOD_H__
#define __IRVING_KIRKWOOD_H__

#include <pyglasstools/Calculator.h> 
#include "cgfunc/CoarseGrainFunction.h"

class PYBIND11_EXPORT IrvingKirkwood : public Calculator
{
    public:
        std::map< std::string, std::shared_ptr< CoarseGrainedField > > m_observables;
        
        IrvingKirkwood( std::shared_ptr< ParticleSystem > sysdata,
                        std::shared_ptr< PairPotential > potential, 
                        std::shared_ptr< CoarseGrainFunction > cgfunc )
            : Calculator(sysdata,potential), m_cgfunc(cgfunc)
        {}; 
        IrvingKirkwood( std::shared_ptr< ParticleSystem > sysdata,
                        std::shared_ptr< PairPotential > potential, 
                        std::shared_ptr< CoarseGrainFunction > cgfunc,
                        std::vector< Eigen::Vector3d > gridpoints )
            : Calculator(sysdata,potential), m_cgfunc(cgfunc), m_gridpoints(gridpoints)
        {
            //We could then distribute the vector here;
        }; 
        ~IrvingKirkwood(){};
        
        void compute();
        
        void setGridpoints(const std::vector< Eigen::Vector3d >& gridpoints)
        {
            m_gridpoints = gridpoints;
        } 

        virtual void addObservable(const std::shared_ptr<CoarseGrainedField>& obs)
        {
            m_observables.insert(std::pair<std::string, std::shared_ptr<CoarseGrainedField> >(obs->name, obs));
        }

        virtual void clearState(unsigned int grid_id)
        {
            for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
                    it->second->clear(grid_id);
        } 

        virtual void computeLocalObsPerGrid(const AboriaParticles::value_type& particle_i, 
                                            double cgval, unsigned int grid_id) 
        {
            for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
            {
                if (it->second->islocal)
                    it->second->accumulate(particle_i,cgval,grid_id);
                else 
                    continue;
            }
        }
        virtual void computePairObsPerGrid( const AboriaParticles::value_type& particle_i, 
                                            const AboriaParticles::value_type& particle_j,
                                            Eigen::Vector3d rij, 
                                            double bondval, unsigned int grid_id)
            {
                for (auto it=m_observables.begin(); it!=m_observables.end(); ++it)
                {
                    if (!it->second->islocal)
                        it->second->accumulate(particle_i,particle_j, rij, m_potential, bondval, grid_id);
                    else 
                        continue;
                }
            }
    private:
        std::shared_ptr< CoarseGrainFunction > m_cgfunc; //!< particle system, equipped with neighbor list
        std::vector< Eigen::Vector3d > m_gridpoints; 
};

//Compute a Global Observable
void IrvingKirkwood::compute()
{
    for (unsigned int i = 0; i < m_gridpoints.size(); ++i)
    {
        clearState(i);
        for( auto p_i = abr::euclidean_search(m_sysdata->particles.get_query(), 
                        abr::vdouble3(m_gridpoints[i][0],m_gridpoints[i][1],m_gridpoints[i][2]), m_cgfunc->getRcut()); 
                        p_i != false; ++p_i)
        {
            //Set grid point X and position of particle ri
            Eigen::Vector3d dr(p_i.dx()[0], p_i.dx()[1],p_i.dx()[2]);
            Eigen::Vector3d x = m_gridpoints[i];
            Eigen::Vector3d ri = m_gridpoints[i]-dr;
            double cgval = m_cgfunc->getDeltaFunc(x,ri);
            computeLocalObsPerGrid(*p_i, cgval, i);
            
            //Compute a list of local observables 
            //Next we loop through the j-th particles for the virial stress
            for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                            abr::get<position>(*p_i), max_rcut); 
                            p_j != false; ++p_j)
            {
                //Make sure the particle is unique
                if (abr::get<abr::id>(*p_i) != abr::get<abr::id>(*p_j))
                {
                    //set the distance between particle i and particle j
                    Eigen::Vector3d rij(-p_j.dx()[0], -p_j.dx()[1], -p_j.dx()[2]);
                    double bondval = m_cgfunc->getBondFunc(x,ri,rij);
                    
                    //Don't forget to set diameters of the potential
                    computePairObsPerGrid(*p_i, *p_j, rij, bondval, i);
                } 
            } 
        }
    }
};


void export_IrvingKirkwood(py::module& m)
{
    py::class_<IrvingKirkwood, Calculator, std::shared_ptr<IrvingKirkwood> >(m,"IrvingKirkwood")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< CoarseGrainFunction > >())
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< CoarseGrainFunction >, std::vector< Eigen::Vector3d > >())
    .def("compute", &IrvingKirkwood::compute)
    .def("setSystemData", &IrvingKirkwood::setSystemData)
    .def("addObservable", &IrvingKirkwood::addObservable)
    ;
};

#endif
