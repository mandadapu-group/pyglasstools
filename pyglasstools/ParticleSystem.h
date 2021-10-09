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
        
        void displaceParticles(); 
        
        void setDisplacement( std::vector< std::vector<double> > atomdisplacement); 
        
        std::vector< std::vector<double> > getDisplacement();
        
        abr::vdouble3 getDisplacementID(unsigned int id);
        
        void setParticleSystemData(  unsigned int numparticles, std::vector<double> atomdiameter, std::vector<double> atommass, 
                                    std::vector< std::vector<double> > atomposition, 
                                    std::vector< std::vector<double> > atomvelocity);

        void setMass(std::vector<double> atommass);
        
        std::vector<double> getMass();
        
        void setDiameter(std::vector<double> atomdiameter);
        
        std::vector<double> getDiameter();
        
        void setAtomPosition(std::vector< std::vector<double> > atomposition);
        
        std::vector< std::vector<double> > getAtomPosition();
        
        void setAtomVelocity(std::vector< std::vector<double> > atomvelocity);
        
        std::vector< std::vector<double> > getAtomVelocity();
        
        std::vector<unsigned int> getNeighbors(std::vector<double> point, double radius);
        
        std::vector< Eigen::Vector3d > getNeighborsDistance(int tag, double radius);
       
        std::shared_ptr<SimBox> simbox;
        AboriaParticles particles;
};

//an export function here
void export_ParticleSystem(py::module& m);
#endif //__SYSTEM_DATA_H__
