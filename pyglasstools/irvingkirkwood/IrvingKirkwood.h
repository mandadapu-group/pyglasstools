#ifndef __IRVING_KIRKWOOD_H__
#define __IRVING_KIRKWOOD_H__

#include <pyglasstools/Calculator.h> 
#include <omp.h>
#include "cgfunc/CoarseGrainFunction.h"

class PYBIND11_EXPORT IrvingKirkwood : public Calculator
{
    
    public:
        IrvingKirkwood( std::shared_ptr< ParticleSystem > sysdata,
                        std::shared_ptr< PairPotential > potential, 
                        std::shared_ptr< CoarseGrainFunction > cgfunc)//,
                        //std::shared_ptr< MPI::Communicator > comm )
            : Calculator(sysdata,potential), m_cgfunc(cgfunc)//, m_comm(comm)
        {
        }; 
        ~IrvingKirkwood(){};
        
        void compute(const std::vector< Eigen::Vector3d >& gridpoints);
        
        virtual void addObservable(const std::shared_ptr<CoarseGrainedField>& obs)
        {
            m_observables.insert(std::pair<std::string, std::shared_ptr<CoarseGrainedField> >(obs->name, obs));
        }
        virtual void printDisplacement()
        {
            for( auto p_i = m_sysdata->particles.begin(); p_i != m_sysdata->particles.end(); ++p_i)
            {
                int id = abr::get<abr::id>(*p_i);
                py::print(abr::get<displacement>(*p_i)[0],abr::get<displacement>(*p_i)[1],id);
                py::print("WHY",abr::get<displacement>(m_sysdata->particles[id])[0],abr::get<displacement>(m_sysdata->particles[id])[1],abr::get<abr::id>(*p_i));
            }
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
        //std::shared_ptr< MPI::Communicator > m_comm;
        std::map< std::string, std::shared_ptr< CoarseGrainedField > > m_observables;
};

//Compute a Global Observable
void IrvingKirkwood::compute(const std::vector< Eigen::Vector3d >& gridpoints)
{
    #pragma omp parallel for
    for (unsigned int i = 0; i < gridpoints.size(); ++i)
    {
        clearState(i);
        for( auto p_i = abr::euclidean_search(m_sysdata->particles.get_query(), 
                        abr::vdouble3(gridpoints[i][0],gridpoints[i][1],gridpoints[i][2]), m_cgfunc->getRcut()); 
                        p_i != false; ++p_i)
        {
            //Set grid point X and position of particle ri
            Eigen::Vector3d dr(p_i.dx()[0], p_i.dx()[1],p_i.dx()[2]);
            Eigen::Vector3d x = gridpoints[i];
            Eigen::Vector3d ri = gridpoints[i]-dr;
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
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< CoarseGrainFunction > >())//, std::shared_ptr< MPI::Communicator >  >())
    .def("compute", &IrvingKirkwood::compute)
    .def("setSystemData", &IrvingKirkwood::setSystemData)
    .def("addObservable", &IrvingKirkwood::addObservable)
    .def("printDisplacement", &IrvingKirkwood::printDisplacement)
    ;
};

#endif
