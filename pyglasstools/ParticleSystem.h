#ifndef __SYSTEM_DATA_H__
#define __SYSTEM_DATA_H__

#include "extern/pybind11/include/pybind11/pybind11.h"
#include <omp.h>
#include "MathAndTypes.h"
#include <pyglasstools/potential/PairPotential.h>
#include "SimBox.h"

#include <Aboria.h>
namespace py = pybind11;
namespace abr = Aboria;

class PYBIND11_EXPORT ParticleSystem
{
    public:
        ParticleSystem( std::shared_ptr< SimBox > simbox, std::shared_ptr<PairPotential> _potential,
                        unsigned int numparticles, std::vector<double> atomdiameter, 
                        std::vector<double> atommass, std::vector< std::vector<double> > atomposition, 
                        std::vector< std::vector<double> > atomvelocity)
                        :   simbox(simbox), potential(_potential), particles(numparticles), 
                            m_numparticles(numparticles)  
        {
            abr::get<diameter>(particles) = atomdiameter;
            abr::get<mass>(particles) = atommass;
            
            #pragma omp parallel for
            for(unsigned int i=0; i < numparticles; i++)
            {
                abr::get<position>(particles[i]) = abr::vdouble3(atomposition[i][0], atomposition[i][1], atomposition[i][2]);
                abr::get<velocity>(particles[i]) = abr::vdouble3(atomvelocity[i][0], atomvelocity[i][1], atomvelocity[i][2]);
            }
            abr::vdouble3 boxmax = abr::vdouble3(simbox->getUpperBound(0), simbox->getUpperBound(1), simbox->getUpperBound(2));
            abr::vdouble3 boxmin = abr::vdouble3(simbox->getLowerBound(0), simbox->getLowerBound(1), simbox->getLowerBound(2));
            abr::vbool3 periodic = abr::vbool3((bool)simbox->getPeriodic(0),(bool)simbox->getPeriodic(1),(bool)simbox->getPeriodic(2));
            particles.init_neighbour_search(boxmin, boxmax, periodic);
        };
        ~ParticleSystem(){};

        void setMass(std::vector<double> atommass)
        {
            if (m_numparticles != atommass.size() )
                throw std::invalid_argument("[ERROR]: Size of mass array mismatch with # of particles!");
            else
                abr::get<mass>(particles) = atommass;
        };
        std::vector<double> getMass()
        {
            return abr::get<mass>(particles); 
        };
        
        void setDiameter(std::vector<double> atomdiameter)
        {
            if (m_numparticles != atomdiameter.size() )
                throw std::invalid_argument("[ERROR]: Size of diameter array mismatch with # of particles!");
            else
                abr::get<diameter>(particles) = atomdiameter;
        };
        std::vector<double> getDiameter()
        {
            return abr::get<diameter>(particles); 
        };
        
        void setAtomPosition(std::vector< std::vector<double> > atomposition)
        {
            if (m_numparticles != atomposition.size() )
                throw std::invalid_argument("[ERROR]: Size of position array mismatch with # of particles!");
            else
            {
                #pragma omp parallel for
                for(unsigned int i=0; i < m_numparticles; i++)
                {
                    abr::get<position>(particles[i]) = abr::vdouble3(atomposition[i][0],atomposition[i][1],atomposition[i][2]);
                }
            }
        };
        void setAtomVelocity(std::vector< std::vector<double> > atomvelocity)
        {
            if (m_numparticles != atomvelocity.size() )
                throw std::invalid_argument("[ERROR]: Size of velocity array mismatch with # of particles!");
            else
            {
                //m_atomvelocity = atomvelocity;
                #pragma omp parallel for
                for(unsigned int i=0; i < m_numparticles; i++)
                {
                    abr::get<velocity>(particles[i]) = abr::vdouble3(atomvelocity[i][0],atomvelocity[i][1],atomvelocity[i][2]);
                }
            }
        };
        
        bool haveVelocities()
        {
            return true; 
        };
        double getScaledRcut()
        {
            //py::print("Rcut from system data", potential->getScaledRcut());
            return potential->getScaledRcut(); 
        };
        
        unsigned int getN()
        {
            return m_numparticles; 
        };
        std::vector<unsigned int> getNeighbors(std::vector<double> point, double radius)
        {

            std::vector<unsigned int> particleID;
            for(    auto particle = abr::euclidean_search(  particles.get_query(), abr::vdouble3(point[0],point[1],point[2]), radius); 
                    particle != false; ++particle)
            {
                particleID.push_back(abr::get<abr::id>(*particle));
            }
            return particleID;
        };
        
        std::shared_ptr<SimBox> simbox;
        std::shared_ptr<PairPotential> potential;
        AboriaParticles particles;
    
    private: 
        unsigned int m_numparticles;
};

//an export function here
void export_ParticleSystem(py::module& m)
{
    py::class_<ParticleSystem, std::shared_ptr<ParticleSystem> >(m,"ParticleSystem")
    .def(py::init<  std::shared_ptr< SimBox >, std::shared_ptr< PairPotential >,
                     unsigned int, std::vector<double>, std::vector<double>, 
                    std::vector< std::vector<double> >, std::vector< std::vector<double> > >())
    .def("getMass", &ParticleSystem::getMass)
    .def("setMass", &ParticleSystem::setMass)
    .def("getDiameter", &ParticleSystem::getDiameter)
    .def("setDiameter", &ParticleSystem::setDiameter)
    .def("setAtomPosition", &ParticleSystem::setAtomPosition)
    .def("setAtomVelocity", &ParticleSystem::setAtomVelocity)
    .def("getNeighbors", &ParticleSystem::getNeighbors)
    .def("getScaledRcut", &ParticleSystem::getScaledRcut)
    ;
};
#endif //__SYSTEM_DATA_H__
