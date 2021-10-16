#ifndef __EIGEN_HESSIAN_H__
#define __EIGEN_HESSIAN_H__

#include "EigenHessianBase.h"

template<int Dim >
class EigenHessian : public EigenHessianBase
{
    public:
        EigenHessian( std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                        std::shared_ptr< EigenManager > manager);
        ~EigenHessian()
        {
            destroyObjects();
        }        
        virtual void destroyObjects()
        {
        };
        
        void assembleObjects();
        
        void setHessianValues(int id_i, int id_j, int real_id, Eigen::Matrix3d offdiag_ij, std::vector< Tripletd > &hessian_triplet);
};

template< int Dim >
EigenHessian<Dim>::EigenHessian(std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                                std::shared_ptr< EigenManager > manager)
    : EigenHessianBase(sysdata, potential, manager)
{
    assembleObjects();
};



//Implementation is still specific to 2D!!!
template< int Dim >
void EigenHessian<Dim>::assembleObjects()
{
        m_manager->notice(5) << "Constructing Hessian Object" << std::endl;
        diagonalsarenonzero = 0;
        
        //Compute the maximum cut-off radius. This is set by either the force cut-off radius or the diameter
        double maxforce_rcut = m_potential->scaled_rcut;
        double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                                std::end(abr::get<diameter>(m_sysdata->particles)) );
        maxforce_rcut *= maxdiameter;
        max_rcut = std::max(maxforce_rcut,maxdiameter); 
        
        //Construct Hessian
        hessian_length = (unsigned int)Dim*m_sysdata->particles.size();
    
        //Construct the hessian matrix, but don't assemble!
        unsigned int estimate_nonzero_entries = (unsigned int)Dim*hessian_length*6;
        std::vector< Tripletd > hessian_triplet;
        hessian_triplet.reserve(estimate_nonzero_entries);
        hessian.resize(hessian_length,hessian_length);
        
        //We loop through the rows owned by the the processor
        std::set< int > particleid_nonneigh;
        
        for (unsigned int i = 0; i < hessian_length; ++i) 
        {
            double neighbors_count = 0;
            //Instantiate the iterator at a particular position
            auto p_i = m_sysdata->particles.begin()+(int)(i/Dim);
            for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                        abr::get<position>(*p_i), max_rcut); p_j != false; ++p_j)
            {
                int id_j = abr::get<abr::id>(*p_j);
                int id_i = abr::get<abr::id>(*p_i);
                
                //We discount processing data for the same particle. 
                if (id_i != id_j)
                {
                    //Set the distance between particle i and particle j. 
                    Eigen::Vector3d rij(p_j.dx()[0], p_j.dx()[1], p_j.dx()[2]);
                    Eigen::Vector3d nij = rij.normalized();
                    
                    double di =  abr::get<diameter>(*p_i);
                    double dj =  abr::get<diameter>(*p_j);
                       
                    //Make sure the particle is unique
                    if (m_potential->getRcut(rij, di, dj) > rij.dot(rij)) 
                    {
                        double factor = m_potential->getBondStiffness(rij, di, dj)+m_potential->getPairForceDivR(rij, di, dj);
                        Eigen::Matrix3d offdiag_ij = -factor*nij*nij.transpose()
                                                     +Eigen::Matrix3d::Identity()*m_potential->getPairForceDivR(rij, di, dj);
                        setHessianValues(id_i, id_j, i, offdiag_ij, hessian_triplet);
                        neighbors_count += 1;
                    }
                }
            }
            if (neighbors_count < 1)
            {
                //Somehow a particle has no neighbors within the interaction cut-off. Add this to the list
                m_manager->notice(7) << "Found particle with no neighbors!" << " \n" << std::string("Particle ID: ")+std::to_string((int)(i/Dim)) << std::endl;
                particleid_nonneigh.insert((int)(i/Dim));
                diagonalsarenonzero = 1;
            }
        }

        hessian.setFromTriplets(hessian_triplet.begin(), hessian_triplet.end());
        if (diagonalsarenonzero > 0)
        {
            m_manager->notice(0) << "[WARNING] Found particles with no neighbors. Any normal mode analysis will be immediately aborted" << std::endl;
            for (auto pid : particleid_nonneigh)
            {
                m_manager->notice(0) << "Particle ID: " << pid << " detected by Process [" << 0 << "]" << std::endl;
            } 
        }
};

template< int Dim >
void EigenHessian<Dim>::setHessianValues(int id_i, int id_j, int real_id, Eigen::Matrix3d offdiag_ij, std::vector< Tripletd > &hessian_triplet)
{
    //Each "set value" must be carefully considered because we don't know the 
    //parallel layout of the matrix beforehand !! But we are allowerd to fill them row-by-row
    //row must also consider indexing at particular 
    //x-component of the row 
    if (Dim*id_i == real_id)
    {
        //The off diagonal term
        hessian_triplet.push_back( Tripletd(Dim*id_i, Dim*id_j, offdiag_ij(0,0)) );
        hessian_triplet.push_back( Tripletd(Dim*id_i, Dim*id_j+1, offdiag_ij(0,1)) );
        //The main diagonal term
        hessian_triplet.push_back( Tripletd(Dim*id_i, Dim*id_i, -offdiag_ij(0,0)) );
        hessian_triplet.push_back( Tripletd(Dim*id_i, Dim*id_i+1, -offdiag_ij(0,1)) );
        if (Dim == 3)
        {
            hessian_triplet.push_back( Tripletd(Dim*id_i, Dim*id_j+2, offdiag_ij(0,2)) );
            hessian_triplet.push_back( Tripletd(Dim*id_i, Dim*id_i+2, -offdiag_ij(0,2)) );
        }
    }
    //y-component of the row
    else if (Dim*id_i+1 == real_id)
    {
        //The off diagonal term
        hessian_triplet.push_back( Tripletd(Dim*id_i+1, Dim*id_j, offdiag_ij(1,0)) );
        hessian_triplet.push_back( Tripletd(Dim*id_i+1, Dim*id_j+1, offdiag_ij(1,1)) );
        //The main diagonal term
        hessian_triplet.push_back( Tripletd(Dim*id_i+1, Dim*id_i, -offdiag_ij(1,0)) );
        hessian_triplet.push_back( Tripletd(Dim*id_i+1, Dim*id_i+1, -offdiag_ij(1,1)) );
        if (Dim == 3)
        {
            hessian_triplet.push_back( Tripletd(Dim*id_i+1, Dim*id_j+2, offdiag_ij(1,2)) );
            hessian_triplet.push_back( Tripletd(Dim*id_i+1, Dim*id_i+2, -offdiag_ij(1,2)) );
        }
    }
    //y-component of the row
    else if (Dim*id_i+2 == real_id && Dim == 3)
    {
        //The off diagonal term
        hessian_triplet.push_back( Tripletd(Dim*id_i+2, Dim*id_j, offdiag_ij(2,0)) );
        hessian_triplet.push_back( Tripletd(Dim*id_i+2, Dim*id_j+1, offdiag_ij(2,1)) );
        hessian_triplet.push_back( Tripletd(Dim*id_i+2, Dim*id_j+2, offdiag_ij(2,2)) );
        
        //the main diagonal term    
        hessian_triplet.push_back( Tripletd(Dim*id_i+2, Dim*id_i, -offdiag_ij(2,0)) );
        hessian_triplet.push_back( Tripletd(Dim*id_i+2, Dim*id_i+1, -offdiag_ij(2,1)) );
        hessian_triplet.push_back( Tripletd(Dim*id_i+2, Dim*id_i+2, -offdiag_ij(2,2)) );
    }
};

template< class T >
void export_EigenHessian(py::module& m, const std::string& name)
{
    py::class_<T, EigenHessianBase, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< EigenManager > >())
    .def("destroyObjects", &T::destroyObjects)
    .def("assembleObjects", &T::assembleObjects)
    .def("setSystemData", &T::setSystemData)
    .def("areDiagonalsNonZero", &T::areDiagonalsNonZero)
    ;
};

#endif
