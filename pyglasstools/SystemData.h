#ifndef __SYSTEM_DATA_H__
#define __SYSTEM_DATA_H__

//#include "extern/aabbcc/src/AABB.h"
#include <Eigen/Dense>
#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/stl.h"
#include <Aboria.h>
namespace py = pybind11;
using namespace Aboria;

class PYBIND11_EXPORT SystemData
{
    public:
        ABORIA_VARIABLE(velocity, vdouble3, "velocity")
        ABORIA_VARIABLE(mass, double, "mass")
        ABORIA_VARIABLE(diameter, double, "diameter")
        typedef Particles< std::tuple<velocity, diameter, mass> >::position position;        
        
        SystemData( unsigned int numparticles, std::vector<double> atomdiameter, 
                    std::vector<double> atommass, std::vector< std::vector<double> > atomposition, 
                    std::vector< std::vector<double> > atomvelocity, std::shared_ptr< SimBox > simbox)
        :   m_numparticles(numparticles), m_particles(numparticles), m_atomposition(atomposition), m_atomvelocity(atomvelocity), m_simbox(simbox) 
        {
            get<diameter>(m_particles) = atomdiameter;
            get<mass>(m_particles) = atommass;
            for(unsigned int i=0; i < numparticles; i++)
            {
                get<position>(m_particles[i]) = vdouble3(atomposition[i][0],atomposition[i][1],atomposition[i][2]);
                get<velocity>(m_particles[i]) = vdouble3(atomvelocity[i][0],atomvelocity[i][1],atomvelocity[i][2]);
            }
            vdouble3 boxmax = vdouble3(m_simbox->getLmax()[0],m_simbox->getLmax()[1],m_simbox->getLmax()[2]);
            vdouble3 boxmin = vdouble3(m_simbox->getLmin()[0],m_simbox->getLmin()[1],m_simbox->getLmin()[2]);
            vbool3 periodic = vbool3(m_simbox->getPeriodic()[0],m_simbox->getPeriodic()[1],m_simbox->getPeriodic()[2]);
            m_particles.init_neighbour_search(boxmin, boxmax, periodic);
        };
        ~SystemData(){};

        void setMass(std::vector<double> atommass)
        {
            if (m_numparticles != atommass.size() )
                throw std::invalid_argument("[ERROR]: Size of mass array mismatch with # of particles!");
            else
                get<mass>(m_particles) = atommass;
        };
        std::vector<double> getMass()
        {
            return get<mass>(m_particles); 
        };
        
        void setDiameter(std::vector<double> atomdiameter)
        {
            if (m_numparticles != atomdiameter.size() )
                throw std::invalid_argument("[ERROR]: Size of diameter array mismatch with # of particles!");
            else
                get<diameter>(m_particles) = atomdiameter;
        };
        std::vector<double> getDiameter()
        {
            return get<diameter>(m_particles); 
        };
        
        void setAtomPosition(std::vector< std::vector<double> > atomposition)
        {
            if (m_numparticles != atomposition.size() )
                throw std::invalid_argument("[ERROR]: Size of position array mismatch with # of particles!");
            else
            {
                m_atomposition = atomposition;
                for(unsigned int i=0; i < m_numparticles; i++)
                {
                    get<position>(m_particles[i]) = vdouble3(atomposition[i][0],atomposition[i][1],atomposition[i][2]);
                }
            }
        };
        
        std::vector< std::vector<double> > getAtomPosition()
        {
            return m_atomposition; 
        };
        
        void setAtomVelocity(std::vector< std::vector<double> > atomvelocity)
        {
            if (m_numparticles != atomvelocity.size() )
                throw std::invalid_argument("[ERROR]: Size of velocity array mismatch with # of particles!");
            else
            {
                m_atomvelocity = atomvelocity;
                for(unsigned int i=0; i < m_numparticles; i++)
                {
                    get<velocity>(m_particles[i]) = vdouble3(atomvelocity[i][0],atomvelocity[i][1],atomvelocity[i][2]);
                }
            }
        };

        std::vector< std::vector<double> > getAtomVelocity()
        {
            return m_atomvelocity; 
        };
        
        unsigned int getN()
        {
            return m_numparticles; 
        };
        std::vector<unsigned int> getNeighbors(std::vector<double> point, double radius)
        {

            std::vector<unsigned int> particleID;
            for (auto particle = euclidean_search(m_particles.get_query(), vdouble3(point[0],point[1],point[2]), radius); particle != false; ++particle)
            {
                particleID.push_back(get<id>(*particle));
            }
            return particleID;
        };
    private: 
        //Atomic properties, fed in to the system
        unsigned int m_numparticles;
        Particles< std::tuple<velocity, diameter, mass> > m_particles;
        
        //positions and velocities in terms of vectors
        std::vector< std::vector<double> > m_atomposition;        
        std::vector< std::vector<double> > m_atomvelocity;        
        //The Simulation Box
        std::shared_ptr<SimBox> m_simbox; 
};

//an export function here
void export_SystemData(py::module& m)
{
    py::class_<SystemData, std::shared_ptr<SystemData> >(m,"SystemData")
    .def(py::init<  unsigned int, std::vector<double>, std::vector<double>, 
                    std::vector< std::vector<double> >, std::vector< std::vector<double> >, std::shared_ptr< SimBox > >())
    .def("getMass", &SystemData::getMass)
    .def("setMass", &SystemData::setMass)
    .def("getDiameter", &SystemData::getDiameter)
    .def("setDiameter", &SystemData::setDiameter)
    .def("getAtomPosition", &SystemData::getAtomPosition)
    .def("setAtomPosition", &SystemData::setAtomPosition)
    .def("getAtomVelocity", &SystemData::getAtomVelocity)
    .def("setAtomVelocity", &SystemData::setAtomVelocity)
    .def("getNeighbors", &SystemData::getNeighbors)
    ;
};
#endif //__SYSTEM_DATA_H__
