#ifndef __SLEPC_HESSIAN_H__
#define __SLEPC_HESSIAN_H__

#include "PETScHessianBase.h"

/*
 * Derived template class for constructing Hessian using PETSc's matrices and arrays.
 * with the specific applications of performing normal mode analysis using SLEPc
 * and computing non-affine elasticity tensors.
 *
 * Template argument Dim is to denote physical dimension of the system, e.g., 2D or 3D. The rationale
 * of making Dim as a template argument is to avoid working with pointers of PETSc vectors and arrays
 * that depend on physical dimensionality of the system.
 */
template<int Dim >
class SLEPcHessian : public PETScHessianBase
{
    public:
        std::string m_hessian_mode;
        SLEPcHessian(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                        std::shared_ptr< PETScManager > manager, std::shared_ptr< MPI::ParallelCommunicator > comm, std::string hessian_mode);
        
        /* Default destructor calls helper function */
        ~SLEPcHessian()
        {
            destroyObjects();
        }        
        /* 
         * Helper function to be used by the class destructor. 
         * It calls in destructor routines from PETSc
         */ 
        void destroyObjects()
        {
            m_manager->printPetscNotice(5,"Destroying PETSc Hessian Object \n");
            MatDestroy(&hessian);
        };
        
        void assembleObjects();
        
        Eigen::VectorXd multiply(const Eigen::VectorXd& vec)
        {
            Vec temp, temp2;
            Eigen::VectorXd out_val = Eigen::VectorXd::Zero(vec.size());
            MatCreateVecs(hessian,NULL,&temp);
            MatCreateVecs(hessian,NULL,&temp2);
            PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(hessian, &Istart, &Iend);
            for(int i = Istart; i < Iend; ++i)
            {
                VecSetValue(temp,i,vec[i], INSERT_VALUES);

            }
            VecAssemblyBegin(temp);
            VecAssemblyEnd(temp);
            VecAssemblyBegin(temp2);
            VecAssemblyEnd(temp2);
            MatMult(hessian,temp,temp2);
            double val;
            for(int i = Istart; i < Iend; ++i)
            {
                VecGetValues(temp2,1,&i, &val);
                out_val[i] = val;
            }
            VecDestroy(&temp);
            VecDestroy(&temp2);
            return out_val;
        }  
        /*
         * Helper function to set the component of the Hessian matrix
         */ 
        void setHessianValues(PetscInt id_i, PetscInt id_j, PetscInt real_id, Eigen::Matrix3d offdiag_ij);
        
        /*
         * Helper function to set the component of the mismatch force vector, denoted as Xi in many papers.
         */ 
        void setMisforceVectorValues(PetscInt id_i,PetscInt real_id,double factor,Eigen::Vector3d rij, Eigen::Vector3d nij);
};


template< int Dim >
SLEPcHessian<Dim>::SLEPcHessian(std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                                std::shared_ptr< PETScManager > manager, std::shared_ptr< MPI::ParallelCommunicator > comm, std::string hessian_mode)
    : PETScHessianBase(sysdata, potential, manager, comm), m_hessian_mode(hessian_mode)
{
    assembleObjects();
};

/*
 * Helper function for assembling the Hessian matrix and the mismatch force vector.
 * This is where both quantities are computed from the particle system data for the first time!
 */    
template< int Dim >
void SLEPcHessian<Dim>::assembleObjects()
{
        //Add command line options 
        ierr = PetscOptionsInsertString(NULL,m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
        m_manager->printPetscNotice(5,"Constructing PETSc Hessian Object \n");
        diagonalsarenonzero = 0;
        
        //Compute the maximum cut-off radius. This is set by either the force cut-off radius or the diameter
        double maxforce_rcut = m_potential->scaled_rcut;
        double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                                std::end(abr::get<diameter>(m_sysdata->particles)) );
        maxforce_rcut *= maxdiameter;
        max_rcut = std::max(maxforce_rcut,maxdiameter); 
        
        hessian_length = (unsigned int)Dim*m_sysdata->particles.size();
        //Construct the PETSc hessian matrix, but don't assemble yet!
        MatCreate(PETSC_COMM_WORLD,&hessian);
        MatSetSizes(hessian,PETSC_DETERMINE,PETSC_DETERMINE,hessian_length,hessian_length);
        MatSetType(hessian,MATAIJ);
        MatSetUp(hessian);
        
        PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(hessian, &Istart, &Iend);
        m_manager->printPetscNotice(5,"Assembling PETSc Sparse Matrix \n");
        
        std::set< int > particleid_nonneigh;
        
        //The first-loop is over the particle index and Cartesian component, as owned by the processor
        for (PetscInt i = Istart; i < Iend; ++i) 
        {
            double neighbors_count = 0;
            
            //Instantiate the iterator at a particular position
            auto p_i = m_sysdata->particles.begin()+(int)(i/Dim); //Divide by i/Dim to get the actual particle index
            
            //The second for-loop is to go through all possible neighbors of the i-th particle
            for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                        abr::get<position>(*p_i), max_rcut); p_j != false; ++p_j)
            {
                PetscInt id_j = abr::get<abr::id>(*p_j);
                PetscInt id_i = abr::get<abr::id>(*p_i);
               
                //We discount processing data for the same particle. 
                if (id_i != id_j)
                {
                    //Set the distance between particle i and particle j. 
                    Eigen::Vector3d rij(p_j.dx()[0], p_j.dx()[1], p_j.dx()[2]);
                    Eigen::Vector3d nij = rij.normalized();
                   
                    //Diameter of each particle 
                    double di =  abr::get<diameter>(*p_i);
                    double dj =  abr::get<diameter>(*p_j);
                       
                    //Make sure the particle is unique
                    if (m_potential->getRcut(rij, di, dj) > rij.dot(rij)) 
                    {
                        //Note PairForceDivR is -phi_r(r)/r, hence the negative sign is already encoded inside factor
                        //double factor = m_potential->getPairForceDivR(rij, di, dj);
                        //double factor = m_potential->getBondStiffness(rij, di, dj);//+m_potential->getPairForceDivR(rij, di, dj);
                        
                        Eigen::Matrix3d offdiag_ij;
                        double factor;
                        //double smallfactor = m_manager->small_factor;
                        if (this->m_hessian_mode == "normal")
                        {
                            //factor = m_potential->getBondStiffness(rij, di, dj)+smallfactor*m_potential->getPairForceDivR(rij, di, dj);
                            //offdiag_ij = -factor*nij*nij.transpose()
                                                         //+smallfactor*Eigen::Matrix3d::Identity()*m_potential->getPairForceDivR(rij, di, dj);
                            factor = m_potential->getBondStiffness(rij, di, dj)+m_potential->getPairForceDivR(rij, di, dj);
                            offdiag_ij = -factor*nij*nij.transpose()
                                                         +Eigen::Matrix3d::Identity()*m_potential->getPairForceDivR(rij, di, dj);
                        }
                        else if (this->m_hessian_mode == "nostress")
                        {
                            factor = m_potential->getBondStiffness(rij, di, dj);
                            offdiag_ij = -factor*nij*nij.transpose();
                        }
                        else if (this->m_hessian_mode == "purestress")
                        {
                            factor = m_potential->getPairForceDivR(rij, di, dj);
                            offdiag_ij = -factor*nij*nij.transpose()+Eigen::Matrix3d::Identity()*factor;
                        }
                        setHessianValues(id_i, id_j, i, offdiag_ij);
                        
                        //setMisforceVectorValues(id_i,i,factor,rij,nij);
                        
                        //We keep increment counter for every neighbor with non-zero interactions
                        neighbors_count += 1;
                    }
                }
            }
            if (neighbors_count < 1)
            {
                std::stringstream string_stream;
                //Somehow a particle has no neighbors within the interaction cut-off. Add this to the list
                string_stream << "Found particle with no neighbors!" << " \n" << std::string("Particle ID: ")+std::to_string((int)(i/Dim)) << std::endl;
                m_manager->printPetscNotice(7,string_stream.str());
                particleid_nonneigh.insert((int)(i/Dim));
                diagonalsarenonzero = 1;
            }
        }
        //Now it's time to gather all the results detected
        std::vector< int > listdiagonalsarenonzero;
        m_comm->all_gather_v(diagonalsarenonzero, listdiagonalsarenonzero);
        diagonalsarenonzero = std::accumulate(listdiagonalsarenonzero.begin(), listdiagonalsarenonzero.end(), 0);
        
        //Assembling all defined PETSc matrices
        MatAssemblyBegin(hessian,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(hessian,MAT_FINAL_ASSEMBLY);
        
        if (diagonalsarenonzero > 0)
        {
            m_manager->printPetscNotice(0,"[WARNING] Found particles with no neighbors. Any normal mode analysis will be immediately aborted \n");
            for (auto pid : particleid_nonneigh)
            {
                m_manager->notice(0) << "Particle ID: " << pid << " detected by Process [" << m_comm->getRank() << "]" << std::endl;
            } 
        }
};

/*
 * Helper function to set the component of the Hessian matrix. Args:
 * id_i: the i-th particle index on the row
 * id_j: the j-th particle index on the column
 * real_id: actual index on the matrix
 * offdiag_ij: the a Dim x Dim matrix representing off-diagonal component of the Hessian matrix. 
 */ 
template< int Dim >
void SLEPcHessian<Dim>::setHessianValues(PetscInt id_i, PetscInt id_j, PetscInt real_id,Eigen::Matrix3d offdiag_ij)
{
    //Each "set value" must be carefully considered because we don't know the 
    //parallel layout of the matrix beforehand !! But we are allowerd to fill them row-by-row
    //row must also consider indexing at particular 
    //x-component of the row 
    if (Dim*id_i == real_id)
    { 
        //The off diagonal term
        MatSetValue(hessian,Dim*id_i, Dim*id_j, offdiag_ij(0,0),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i, Dim*id_j+1, offdiag_ij(0,1),ADD_VALUES);
        
        //The main diagonal term
        MatSetValue(hessian,Dim*id_i, Dim*id_i, -offdiag_ij(0,0),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i, Dim*id_i+1, -offdiag_ij(0,1),ADD_VALUES);
         
        if (Dim == 3)
        {
            //The off diagonal term
            MatSetValue(hessian,Dim*id_i, Dim*id_j+2, offdiag_ij(0,2),ADD_VALUES);
            
            //The main diagonal term
            MatSetValue(hessian,Dim*id_i, Dim*id_i+2, -offdiag_ij(0,2),ADD_VALUES);
        }
    }
    //y-component of the row
    else if (Dim*id_i+1 == real_id)
    {
        //The off diagonal term
        MatSetValue(hessian,Dim*id_i+1, Dim*id_j, offdiag_ij(1,0),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i+1, Dim*id_j+1, offdiag_ij(1,1),ADD_VALUES);
        
        //the main diagonal term    
        MatSetValue(hessian,Dim*id_i+1, Dim*id_i, -offdiag_ij(1,0),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i+1, Dim*id_i+1, -offdiag_ij(1,1),ADD_VALUES);
        if (Dim == 3)
        {
            //The off diagonal term
            MatSetValue(hessian,Dim*id_i+1, Dim*id_j+2, offdiag_ij(1,2),ADD_VALUES);
            
            //the main diagonal term    
            MatSetValue(hessian,Dim*id_i+1, Dim*id_i+2, -offdiag_ij(1,2),ADD_VALUES);
        }
    }
    //y-component of the row
    else if (Dim*id_i+2 == real_id && Dim == 3)
    {
        //The off diagonal term
        MatSetValue(hessian,Dim*id_i+2, Dim*id_j, offdiag_ij(2,0),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i+2, Dim*id_j+1, offdiag_ij(2,1),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i+2, Dim*id_j+2, offdiag_ij(2,2),ADD_VALUES);
        
        //the main diagonal term    
        MatSetValue(hessian,Dim*id_i+2, Dim*id_i, -offdiag_ij(2,0),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i+2, Dim*id_i+1, -offdiag_ij(2,1),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i+2, Dim*id_i+2, -offdiag_ij(2,2),ADD_VALUES);
    }
};

/*
 * Helper function to export SLEPcHessian to Python
 */
template< class T >
void export_SLEPcHessian(py::module& m, const std::string& name)
{
    py::class_<T, PETScHessianBase, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< PETScManager > , std::shared_ptr< MPI::ParallelCommunicator >, std::string >())
    .def("destroyObjects", &T::destroyObjects)
    .def("assembleObjects", &T::assembleObjects)
    .def("setSystemData", &T::setSystemData)
    .def("multiply",&T::multiply)
    .def("areDiagonalsNonZero", &T::areDiagonalsNonZero)
    .def_readwrite("m_hessian_mode",&T::m_hessian_mode)
    ;
};

#endif
