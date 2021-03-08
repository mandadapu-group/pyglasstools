#ifndef __SYSTEM_DATA_H__
#define __SYSTEM_DATA_H__

#include <pybind11/pybind11.h>
#include <omp.h>
#include "MathAndTypes.h"
#include "SimBox.h"
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <Aboria.h>
namespace py = pybind11;
namespace abr = Aboria;

class PYBIND11_EXPORT ParticleSystem
{
    public:
        ParticleSystem( std::shared_ptr< SimBox > simbox,
                        unsigned int numparticles, std::vector<double> atomdiameter, 
                        std::vector<double> atommass, std::vector< std::vector<double> > atomposition, 
                        std::vector< std::vector<double> > atomvelocity);
        ~ParticleSystem(){};
        
        void updateParticleSystem(  unsigned int numparticles, std::vector<double> atomdiameter, std::vector<double> atommass, 
                                    std::vector< std::vector<double> > atomposition, 
                                    std::vector< std::vector<double> > atomvelocity);
        
        void moveParticles( std::vector< std::vector<double> > atomposition); 
        
        void moveParticles(); 
        
        void setDisplacement( std::vector< std::vector<double> > atomdisplacement); 
        
        std::vector< std::vector<double> > getDisplacement();
        
        abr::vdouble3 getDisplacementID(unsigned int id);
        
        void setParticleSystemData(  unsigned int numparticles, std::vector<double> atomdiameter, std::vector<double> atommass, 
                                    std::vector< std::vector<double> > atomposition, 
                                    std::vector< std::vector<double> > atomvelocity);

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
        
        std::vector< Eigen::Vector3d > getNeighborsDistance(int tag, double radius)
        {

            std::vector< Eigen::Vector3d > rij_list;
            auto p_i = particles.begin()+tag;
            for( auto p_j = abr::euclidean_search(particles.get_query(), 
                        abr::get<position>(*p_i), radius); p_j != false; ++p_j)
            {
                Eigen::Vector3d dr(-p_j.dx()[0], -p_j.dx()[1],-p_j.dx()[2]);
                rij_list.push_back(dr);
            }
            return rij_list;
        };
        
        std::shared_ptr<SimBox> simbox;
        AboriaParticles particles;
};

//an export function here
void export_ParticleSystem(py::module& m);
#endif //__SYSTEM_DATA_H__
