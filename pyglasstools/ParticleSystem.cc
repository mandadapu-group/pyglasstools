#include "ParticleSystem.h"

ParticleSystem::ParticleSystem( std::shared_ptr< SimBox > simbox,
                                unsigned int numparticles, std::vector<double> atomdiameter, 
                                std::vector<double> atommass, std::vector< std::vector<double> > atomposition, 
                                std::vector< std::vector<double> > atomvelocity)
                                :   simbox(simbox), particles(numparticles) 
{
    setParticleSystemData(numparticles, atomdiameter, atommass, atomposition, atomvelocity);
    abr::vdouble3 boxmax = abr::vdouble3(simbox->upperbound[0], simbox->upperbound[1], simbox->upperbound[2]);
    abr::vdouble3 boxmin = abr::vdouble3(simbox->lowerbound[0], simbox->lowerbound[1], simbox->lowerbound[2]);
    abr::vbool3 periodic = abr::vbool3((bool)simbox->periodic[0],(bool)simbox->periodic[1],(bool)simbox->periodic[2]);
    particles.init_neighbour_search(boxmin, boxmax, periodic);
};

void ParticleSystem::setParticleSystemData(  unsigned int numparticles, std::vector<double> atomdiameter, std::vector<double> atommass, 
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

void ParticleSystem::updateParticleSystem(  unsigned int numparticles, std::vector<double> atomdiameter, std::vector<double> atommass, 
                            std::vector< std::vector<double> > atomposition, 
                            std::vector< std::vector<double> > atomvelocity)
{
    setParticleSystemData(numparticles, atomdiameter, atommass, atomposition, atomvelocity);
    particles.update_positions();
};

void ParticleSystem::moveParticles( std::vector< std::vector<double> > atomposition) 
{
    #pragma omp parallel for
    for(unsigned int i=0; i < particles.size(); i++)
    {
        abr::get<position>(particles[i]) += abr::vdouble3(atomposition[i][0], atomposition[i][1], atomposition[i][2]);
    }
    particles.update_positions();
};

void ParticleSystem::displaceParticles() 
{
    #pragma omp parallel for
    for(unsigned int i=0; i < particles.size(); i++)
    {
        abr::get<position>(particles[i]) += abr::get<displacement>(particles[i]);
    }
    particles.update_positions();
};

void ParticleSystem::setDisplacement( std::vector< std::vector<double> > atomdisplacement) 
{
    #pragma omp parallel for
    for(unsigned int i=0; i < particles.size(); i++)
    {
        abr::get<displacement>(particles[i]) = abr::vdouble3(atomdisplacement[i][0], atomdisplacement[i][1], atomdisplacement[i][2]);
    }
};

std::vector< std::vector<double> > ParticleSystem::getDisplacement()
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

abr::vdouble3 ParticleSystem::getDisplacementID(unsigned int id)
{
    return abr::get<displacement>(particles[id]);
};


void ParticleSystem::setMass(std::vector<double> atommass)
{
    if (particles.size() != atommass.size() )
        throw std::invalid_argument("[ERROR]: Size of mass array mismatch with # of particles!");
    else
        abr::get<mass>(particles) = atommass;
};

std::vector<double> ParticleSystem::getMass()
{
    return abr::get<mass>(particles); 
};

void ParticleSystem::setDiameter(std::vector<double> atomdiameter)
{
    if (particles.size() != atomdiameter.size() )
        throw std::invalid_argument("[ERROR]: Size of diameter array mismatch with # of particles!");
    else
        abr::get<diameter>(particles) = atomdiameter;
};

std::vector<double> ParticleSystem::getDiameter()
{
    return abr::get<diameter>(particles); 
};

void ParticleSystem::setAtomPosition(std::vector< std::vector<double> > atomposition)
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

std::vector< std::vector<double> > ParticleSystem::getAtomPosition()
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


void ParticleSystem::setAtomVelocity(std::vector< std::vector<double> > atomvelocity)
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

std::vector< std::vector<double> > ParticleSystem::getAtomVelocity()
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

std::vector<unsigned int> ParticleSystem::getNeighbors(std::vector<double> point, double radius)
{

    std::vector<unsigned int> particleID;
    for(    auto particle = abr::euclidean_search(  particles.get_query(), abr::vdouble3(point[0],point[1],point[2]), radius); 
            particle != false; ++particle)
    {
        particleID.push_back(abr::get<abr::id>(*particle));
    }
    return particleID;
};

std::vector< Eigen::Vector3d > ParticleSystem::getNeighborsDistance(int tag, double radius)
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
    .def("getNeighborsDistance", &ParticleSystem::getNeighborsDistance)
    .def("update", &ParticleSystem::updateParticleSystem) 
    ;
};
