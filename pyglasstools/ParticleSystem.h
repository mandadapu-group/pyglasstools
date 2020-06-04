#ifndef __SYSTEM_DATA_H__
#define __SYSTEM_DATA_H__

#include "extern/pybind11/include/pybind11/pybind11.h"
#include <omp.h>
#include "MathAndTypes.h"
#include "SimBox.h"

#include <Aboria.h>
namespace py = pybind11;
namespace abr = Aboria;

class PYBIND11_EXPORT ParticleSystem
{
    public:
        ParticleSystem( std::shared_ptr< SimBox > simbox,
                        unsigned int numparticles, std::vector<double> atomdiameter, 
                        std::vector<double> atommass, std::vector< std::vector<double> > atomposition, 
                        std::vector< std::vector<double> > atomvelocity)
                        :   simbox(simbox), particles(numparticles) 
        {
            setParticleSystemData(numparticles, atomdiameter, atommass, atomposition, atomvelocity);
            abr::vdouble3 boxmax = abr::vdouble3(simbox->getUpperBound(0), simbox->getUpperBound(1), simbox->getUpperBound(2));
            abr::vdouble3 boxmin = abr::vdouble3(simbox->getLowerBound(0), simbox->getLowerBound(1), simbox->getLowerBound(2));
            abr::vbool3 periodic = abr::vbool3((bool)simbox->getPeriodic(0),(bool)simbox->getPeriodic(1),(bool)simbox->getPeriodic(2));
            particles.init_neighbour_search(boxmin, boxmax, periodic);
        };
        ~ParticleSystem(){};
        void updateParticleSystem(  unsigned int numparticles, std::vector<double> atomdiameter, std::vector<double> atommass, 
                                    std::vector< std::vector<double> > atomposition, 
                                    std::vector< std::vector<double> > atomvelocity)
        {
            setParticleSystemData(numparticles, atomdiameter, atommass, atomposition, atomvelocity);
            particles.update_positions();
        }
        
        void moveParticles( std::vector< std::vector<double> > atomposition) 
        {
            #pragma omp parallel for
            for(unsigned int i=0; i < particles.size(); i++)
            {
                abr::get<position>(particles[i]) += abr::vdouble3(atomposition[i][0], atomposition[i][1], atomposition[i][2]);
            }
            particles.update_positions();
        };
        
        void moveParticles() 
        {
            #pragma omp parallel for
            for(unsigned int i=0; i < particles.size(); i++)
            {
                abr::get<position>(particles[i]) += abr::get<displacement>(particles[i]);
            }
            particles.update_positions();
        };
        
        void setDisplacement( std::vector< std::vector<double> > atomdisplacement) 
        {
            #pragma omp parallel for
            for(unsigned int i=0; i < particles.size(); i++)
            {
                abr::get<displacement>(particles[i]) = abr::vdouble3(atomdisplacement[i][0], atomdisplacement[i][1], atomdisplacement[i][2]);
            }
        };
        
        std::vector< std::vector<double> > getDisplacement()
        {
            std::vector< std::vector<double> > atomdisplacement(particles.size(), std::vector<double>(3,0));
            #pragma omp parallel for
            for(unsigned int i=0; i < particles.size(); i++)
            {
                atomdisplacement[i][0] = abr::get<displacement>(particles[i])[0];
                atomdisplacement[i][1] = abr::get<displacement>(particles[i])[1];
                atomdisplacement[i][2] = abr::get<displacement>(particles[i])[2];
            }
            return atomdisplacement;
        };
        abr::vdouble3 getDisplacementID(unsigned int id)
        {
            return abr::get<displacement>(particles[id]);
        };
        
        void setParticleSystemData(  unsigned int numparticles, std::vector<double> atomdiameter, std::vector<double> atommass, 
                                    std::vector< std::vector<double> > atomposition, 
                                    std::vector< std::vector<double> > atomvelocity)
        {
            abr::get<diameter>(particles) = atomdiameter;
            abr::get<mass>(particles) = atommass;
            
            #pragma omp parallel for
            for(unsigned int i=0; i < numparticles; i++)
            {
                abr::get<position>(particles[i]) = abr::vdouble3(atomposition[i][0], atomposition[i][1], atomposition[i][2]);
                abr::get<velocity>(particles[i]) = abr::vdouble3(atomvelocity[i][0], atomvelocity[i][1], atomvelocity[i][2]);
                abr::get<displacement>(particles[i]) = abr::vdouble3(0,0,0);
            }
        };

        void setMass(std::vector<double> atommass)
        {
            if (particles.size() != atommass.size() )
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
            if (particles.size() != atomdiameter.size() )
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
            if (particles.size() != atomposition.size() )
                throw std::invalid_argument("[ERROR]: Size of position array mismatch with # of particles!");
            else
            {
                #pragma omp parallel for
                for(unsigned int i=0; i < particles.size(); i++)
                {
                    abr::get<position>(particles[i]) = abr::vdouble3(atomposition[i][0],atomposition[i][1],atomposition[i][2]);
                }
            }
        };
        std::vector< std::vector<double> > getAtomPosition()
        {
            std::vector< std::vector<double> > atomposition(particles.size(), std::vector<double>(3,0));
            #pragma omp parallel for
            for(unsigned int i=0; i < particles.size(); i++)
            {
                atomposition[i][0] = abr::get<position>(particles[i])[0];
                atomposition[i][1] = abr::get<position>(particles[i])[1];
                atomposition[i][2] = abr::get<position>(particles[i])[2];
            }
            return atomposition;
        };
        
        void setAtomVelocity(std::vector< std::vector<double> > atomvelocity)
        {
            if (particles.size() != atomvelocity.size() )
                throw std::invalid_argument("[ERROR]: Size of velocity array mismatch with # of particles!");
            else
            {
                //m_atomvelocity = atomvelocity;
                #pragma omp parallel for
                for(unsigned int i=0; i < particles.size(); i++)
                {
                    abr::get<velocity>(particles[i]) = abr::vdouble3(atomvelocity[i][0],atomvelocity[i][1],atomvelocity[i][2]);
                }
                particles.update_positions();
            }
        };
        
        std::vector< std::vector<double> > getAtomVelocity()
        {
            std::vector< std::vector<double> > atomvelocity(particles.size(), std::vector<double>(3,0));
            #pragma omp parallel for
            for(unsigned int i=0; i < particles.size(); i++)
            {
                atomvelocity[i][0] = abr::get<velocity>(particles[i])[0];
                atomvelocity[i][1] = abr::get<velocity>(particles[i])[1];
                atomvelocity[i][2] = abr::get<velocity>(particles[i])[2];
            }
            return atomvelocity;
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
        AboriaParticles particles;
};

//an export function here
void export_ParticleSystem(py::module& m)
{
    py::class_<ParticleSystem, std::shared_ptr<ParticleSystem> >(m,"ParticleSystem")
    .def(py::init<  std::shared_ptr< SimBox >,
                     unsigned int, std::vector<double>, std::vector<double>, 
                    std::vector< std::vector<double> >, std::vector< std::vector<double> > >())
    .def("getMass", &ParticleSystem::getMass)
    .def("setMass", &ParticleSystem::setMass)
    .def("getDiameter", &ParticleSystem::getDiameter)
    .def("setDiameter", &ParticleSystem::setDiameter)
    .def("getDisplacement", &ParticleSystem::getDisplacement)
    .def("setDisplacement", &ParticleSystem::setDisplacement)
    .def("setAtomPosition", &ParticleSystem::setAtomPosition)
    .def("setAtomVelocity", &ParticleSystem::setAtomVelocity)
    .def("getAtomPosition", &ParticleSystem::getAtomPosition)
    .def("getAtomVelocity", &ParticleSystem::getAtomVelocity)
    .def("moveParticles", (void (ParticleSystem::*)(std::vector< std::vector<double> >)) &ParticleSystem::moveParticles)
    .def("moveParticles", (void (ParticleSystem::*)( )) &ParticleSystem::moveParticles)
    .def("getNeighbors", &ParticleSystem::getNeighbors)
    .def("update", &ParticleSystem::updateParticleSystem) 
    ;
};
#endif //__SYSTEM_DATA_H__
