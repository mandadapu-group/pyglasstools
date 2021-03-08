#ifndef __SLEPC_HESSIAN_H__
#define __SLEPC_HESSIAN_H__

#include "HessianBase.h"

template<int Dim >
class SLEPcHessian : public PETScHessianBase
{
    public:
        MatNullSpace constant;
        Vec nullvec[Dim];

        SLEPcHessian(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                        std::shared_ptr< PETScManager > manager, std::shared_ptr< MPI::Communicator > comm);
        ~SLEPcHessian()
        {
            destroyObjects();
        }        
        
        void destroyObjects()
        {
            MatNullSpaceDestroy(&constant);
            for(int i = 0; i < Dim; ++i)
            {
                VecDestroy(&nullvec[i]);
            }
            MatDestroy(&hessian);
            MatDestroy(&misforce);
        };
        
        void assembleObjects();
        
        void setHessianValues(PetscInt id_i, PetscInt id_j, PetscInt Istart, PetscInt Iend, PetscInt real_id,Eigen::Matrix3d offdiag_ij);
        void setMisforceVectorValues(PetscInt id_i,PetscInt real_id,double factor,Eigen::Vector3d rij, Eigen::Vector3d nij);
        void setNullSpaceBasis(PetscInt id_i, PetscInt Istart, PetscInt Iend, PetscInt real_id);
};
//dynamic_pointer_cast<Base, Derived>
template< int Dim >
SLEPcHessian<Dim>::SLEPcHessian(std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                                std::shared_ptr< PETScManager > manager, std::shared_ptr< MPI::Communicator > comm)
    : PETScHessianBase(sysdata, potential, manager, comm)
{
    assembleObjects();
};

//Implementation is still specific to 2D!!!
template< int Dim >
void SLEPcHessian<Dim>::assembleObjects()
{
        //Add command line options
        ierr = PetscOptionsInsertString(NULL,m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
        m_manager->printPetscNotice(5,"Constructing PETSc Hessian Object \n");
        diagonalsarenonzero = 0;
        //Compute the max_rcut 
        double maxforce_rcut = m_potential->scaled_rcut;
        double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                                std::end(abr::get<diameter>(m_sysdata->particles)) );
        maxforce_rcut *= maxdiameter;
        max_rcut = std::max(maxforce_rcut,maxdiameter); //whichever is maximum
        
        //Construct Hessian
        hessian_length = (unsigned int)Dim*m_sysdata->particles.size();
        
        //Construct the hessian matrix, but don't assemble!
        MatCreate(PETSC_COMM_WORLD,&hessian);
        MatSetSizes(hessian,PETSC_DETERMINE,PETSC_DETERMINE,hessian_length,hessian_length);
        MatSetType(hessian,MATAIJ);
        MatSetUp(hessian);
        //set the nullspace
        m_manager->printPetscNotice(5,"Assemble PETSc Sparse Matrix \n");
        
        MatCreate(PETSC_COMM_WORLD,&misforce);
        MatSetSizes(misforce,PETSC_DETERMINE,PETSC_DETERMINE,hessian_length,(unsigned int)Dim*((unsigned int)Dim+1)/2);
        MatSetType(misforce,MATDENSE);
        MatSetUp(misforce);
        
        PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(hessian, &Istart, &Iend);
        
        //We loop through the rows owned by the the processor
        std::set< int > particleid_nonneigh;
        
        for (PetscInt i = Istart; i < Iend; ++i) 
        {
            double neighbors_count = 0;
            //Instantiate the iterator at a particular position
            auto p_i = m_sysdata->particles.begin()+(int)(i/Dim);
            for( auto p_j = abr::euclidean_search(m_sysdata->particles.get_query(), 
                        abr::get<position>(*p_i), max_rcut); p_j != false; ++p_j)
            {
                PetscInt id_j = abr::get<abr::id>(*p_j);
                PetscInt id_i = abr::get<abr::id>(*p_i);
                
                if (id_i != id_j)
                {
                    //Compute a list of local obsercavles 
                    //Next we loop through the j-th particles for the virial stress
                    //set the distance between particle i and particle j
                    Eigen::Vector3d rij(-p_j.dx()[0], -p_j.dx()[1], -p_j.dx()[2]);
                    Eigen::Vector3d nij = rij.normalized();
                    
                    double di =  abr::get<diameter>(*p_i);
                    double dj =  abr::get<diameter>(*p_j);
                       
                    //Make sure the particle is unique
                    if (m_potential->getRcut(rij, di, dj) > rij.dot(rij)) 
                    {
                        //Note PairForce is - Phi_r(r), hence the negative sign is encoded inside factor
                        double factor = m_potential->getBondStiffness(rij, di, dj)+m_potential->getPairForceDivR(rij, di, dj);
                        Eigen::Matrix3d offdiag_ij = -factor*nij*nij.transpose()
                                                     +Eigen::Matrix3d::Identity()*m_potential->getPairForceDivR(rij, di, dj);
                        setHessianValues(id_i, id_j, Istart, Iend, i, offdiag_ij);
                        setMisforceVectorValues(id_i,i,factor,rij,nij);
                        neighbors_count += 1;
                    }
                }
            }
            if (neighbors_count < 1)
            {
                //Somehow a particle has no neighbors within the interaction cut-off. Add this to the list
                m_manager->notice(7,m_comm->getRank()) << "Found particle with no neighbors!" << " \n" << std::string("Particle ID: ")+std::to_string((int)(i/Dim)) << std::endl;
                particleid_nonneigh.insert((int)(i/Dim));
                diagonalsarenonzero = 1;
            }
        }

        //Now it's time to gather all the results detected
        std::vector< int > listdiagonalsarenonzero;
        m_comm->all_gather_v(diagonalsarenonzero, listdiagonalsarenonzero);
        diagonalsarenonzero = std::accumulate(listdiagonalsarenonzero.begin(), listdiagonalsarenonzero.end(), 0);
        
        //std::vector< std::vector< double > > totalparticleid_nonneigh;        
        if (diagonalsarenonzero < 1)
        {
            m_manager->printPetscNotice(5,"Begin Assembling the null space of PETSc matrix\n");
            
            for (int i = 0; i < Dim; ++i)
            {
                MatCreateVecs(hessian,NULL,&nullvec[i]);
            }
            
            for (PetscInt i = Istart; i < Iend; ++i) 
            {
                //Instantiate the iterator at a particular position
                auto p_i = m_sysdata->particles.begin()+(int)(i/Dim);
                PetscInt id_i = abr::get<abr::id>(*p_i);
                setNullSpaceBasis(id_i, Istart, Iend, i);
            }
            
            m_manager->printPetscNotice(5,"Assemble the null space of PETSc matrix\n");
            for (int i = 0; i < Dim; ++i)
            {
                VecAssemblyBegin(nullvec[i]);
                VecAssemblyEnd(nullvec[i]);
                VecNormalize(nullvec[i],NULL);
            }

            MatAssemblyBegin(hessian,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(hessian,MAT_FINAL_ASSEMBLY);
            
            MatAssemblyBegin(misforce,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(misforce,MAT_FINAL_ASSEMBLY);
            
            MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_FALSE,Dim,nullvec,&constant);
            MatSetNullSpace(hessian,constant);
        }
        else
        {
            MatAssemblyBegin(hessian,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(hessian,MAT_FINAL_ASSEMBLY);
            
            MatAssemblyBegin(misforce,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(misforce,MAT_FINAL_ASSEMBLY);
            
            m_manager->printPetscNotice(0,"[WARNING] Found particles with no neighbors. Any normal mode analysis will be immediately aborted \n");
            for (auto pid : particleid_nonneigh)
            {
                m_manager->notice(0) << "Particle ID: " << pid << " detected by Process [" << m_comm->getRank() << "]" << std::endl;
            } 
        }
};

template< int Dim >
void SLEPcHessian<Dim>::setHessianValues(PetscInt id_i, PetscInt id_j, PetscInt Istart, PetscInt Iend, PetscInt real_id,Eigen::Matrix3d offdiag_ij)
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
            MatSetValue(hessian,Dim*id_i, Dim*id_j+2, offdiag_ij(0,2),ADD_VALUES);
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
            MatSetValue(hessian,Dim*id_i+1, Dim*id_j+2, offdiag_ij(1,2),ADD_VALUES);
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

template< int Dim >
void SLEPcHessian<Dim>::setMisforceVectorValues(PetscInt id_i, int real_id, double factor, Eigen::Vector3d rij, Eigen::Vector3d nij)
{
    if (Dim*id_i == real_id)
    { 
        MatSetValue(misforce,Dim*id_i,0, -factor*rij[0]*nij[0]*nij[0],ADD_VALUES); 
        MatSetValue(misforce,Dim*id_i,1, -factor*rij[0]*nij[1]*nij[0],ADD_VALUES); 
        MatSetValue(misforce,Dim*id_i,2, -factor*rij[1]*nij[1]*nij[0],ADD_VALUES); 

        if (Dim == 3)
        {
            MatSetValue(misforce,Dim*id_i,3, -factor*rij[1]*nij[2]*nij[0],ADD_VALUES); 
            MatSetValue(misforce,Dim*id_i,4, -factor*rij[2]*nij[2]*nij[0],ADD_VALUES); 
            MatSetValue(misforce,Dim*id_i,5, -factor*rij[0]*nij[2]*nij[0],ADD_VALUES); 
        }
    }
    
    //y-component of the row
    else if (Dim*id_i+1 == real_id)
    {
        MatSetValue(misforce,Dim*id_i+1,0, -factor*rij[0]*nij[0]*nij[1],ADD_VALUES); 
        MatSetValue(misforce,Dim*id_i+1,1, -factor*rij[0]*nij[1]*nij[1],ADD_VALUES); 
        MatSetValue(misforce,Dim*id_i+1,2, -factor*rij[1]*nij[1]*nij[1],ADD_VALUES); 
        if (Dim == 3)
        {
            MatSetValue(misforce,Dim*id_i+1,3, -factor*rij[1]*nij[2]*nij[1],ADD_VALUES); 
            MatSetValue(misforce,Dim*id_i+1,4, -factor*rij[2]*nij[2]*nij[1],ADD_VALUES); 
            MatSetValue(misforce,Dim*id_i+1,5, -factor*rij[0]*nij[2]*nij[1],ADD_VALUES); 
        }
    }
    //y-component of the row
    else if (Dim*id_i+2 == real_id && Dim == 3)
    {
        MatSetValue(misforce,Dim*id_i+2,0, -factor*rij[0]*nij[0]*nij[2],ADD_VALUES); 
        MatSetValue(misforce,Dim*id_i+2,1, -factor*rij[0]*nij[1]*nij[2],ADD_VALUES); 
        MatSetValue(misforce,Dim*id_i+2,2, -factor*rij[1]*nij[1]*nij[2],ADD_VALUES); 
        MatSetValue(misforce,Dim*id_i+2,3, -factor*rij[1]*nij[2]*nij[2],ADD_VALUES); 
        MatSetValue(misforce,Dim*id_i+2,4, -factor*rij[2]*nij[2]*nij[2],ADD_VALUES); 
        MatSetValue(misforce,Dim*id_i+2,5, -factor*rij[0]*nij[2]*nij[2],ADD_VALUES); 
    }
    
}

template< int Dim >
void SLEPcHessian<Dim>::setNullSpaceBasis(PetscInt id_i, PetscInt Istart, PetscInt Iend, PetscInt real_id)
{
    if (Dim*id_i == real_id)
    { 
        VecSetValue(nullvec[0],Dim*id_i,1, INSERT_VALUES);
        VecSetValue(nullvec[1],Dim*id_i,0, INSERT_VALUES);
        if (Dim == 3)
            VecSetValue(nullvec[Dim-1],Dim*id_i,0, INSERT_VALUES);
    }
    //y-component of the row
    else if (Dim*id_i+1 == real_id)
    {
        VecSetValue(nullvec[0],Dim*id_i+1,0, INSERT_VALUES);
        VecSetValue(nullvec[1],Dim*id_i+1,1, INSERT_VALUES);
        if (Dim == 3)
            VecSetValue(nullvec[Dim-1],Dim*id_i+1,0, INSERT_VALUES);
    }
    
    //y-component of the row
    else if (Dim*id_i+2 == real_id && Dim == 3)
    {
        VecSetValue(nullvec[0],Dim*id_i+2,0, INSERT_VALUES);
        VecSetValue(nullvec[1],Dim*id_i+2,0, INSERT_VALUES);
        if (Dim == 3)
            VecSetValue(nullvec[Dim-1],Dim*id_i+2,1, INSERT_VALUES);
    }
};

template< class T >
void export_SLEPcHessian(py::module& m, const std::string& name)
{
    py::class_<T, PETScHessianBase, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< PETScManager > , std::shared_ptr< MPI::Communicator > >())
    .def("destroyObjects", &T::destroyObjects)
    .def("assembleObjects", &T::assembleObjects)
    .def("setSystemData", &T::setSystemData)
    .def("areDiagonalsNonZero", &T::areDiagonalsNonZero)
    ;
};

#endif
