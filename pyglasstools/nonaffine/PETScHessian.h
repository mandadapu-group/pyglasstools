#ifndef __PETSC_HESSIAN_H__
#define __PETSC_HESSIAN_H__

#include "HessianBase.h"

template< int Dim >
class PETScHessian : public PETScHessianBase
{
    public:
        MatNullSpace constant;
        Vec nullvec[Dim];//, nully
        
        PETScHessian(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                        std::shared_ptr< PETScManager > manager, std::shared_ptr< MPI::Communicator > comm);
        ~PETScHessian()
        {
            destroyPETScObjects();
        }        
        void destroyPETScObjects()
        {
            MatNullSpaceDestroy(&constant);
            for(int i = 0; i < Dim; ++i)
            {
                VecDestroy(&nullvec[i]);
            }
            MatDestroy(&hessian);
            MatDestroy(&misforce);
        };
        void assemblePETScObjects();
        void setHessianValues(PetscInt id_i, PetscInt id_j, PetscInt Istart, PetscInt Iend, PetscInt real_id,Eigen::Matrix3d offdiag_ij);
        void setNullSpaceBasis(PetscInt id_i, PetscInt Istart, PetscInt Iend, PetscInt real_id);
};

//dynamic_pointer_cast<Base, Derived>
template< int Dim >
PETScHessian<Dim>::PETScHessian(std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                                std::shared_ptr< PETScManager > manager, std::shared_ptr< MPI::Communicator > comm)
    : PETScHessianBase(sysdata, potential, manager, comm)
{
    assemblePETScObjects();
};

//Implementation is still specific to 2D!!!
template< int Dim >
void PETScHessian<Dim>::assemblePETScObjects()
{
        //Add command line options
        ierr = PetscOptionsInsertString(NULL,m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
        m_manager->printPetscNotice(5,"Constructing PETSc Hessian Object \n");
       
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
        
        //Construct a "mock" misforce matrix, because we don't really need it
        MatCreate(PETSC_COMM_WORLD,&misforce);
        MatSetSizes(misforce,PETSC_DETERMINE,PETSC_DETERMINE,m_comm->getSizeGlobal(),m_comm->getSizeGlobal());
        MatSetType(misforce,MATDENSE);
        MatSetUp(misforce);
        //set the nullspace
        m_manager->printPetscNotice(5,"Assemble PETSc Sparse Matrix \n");
        
        PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(hessian, &Istart, &Iend);
        //We loop through the rows owned by the the processor
        for (PetscInt i = Istart; i < Iend; ++i) 
        {
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
                        double factor = m_potential->getBondStiffness(rij, di, dj)+m_potential->getPairForce(rij, di, dj);
                        Eigen::Matrix3d offdiag_ij = -factor*nij*nij.transpose()
                                                     +Eigen::Matrix3d::Identity()*m_potential->getPairForce(rij, di, dj);
                        setHessianValues(id_i, id_j, Istart, Iend, i, offdiag_ij);
                    }
                }
            }
        }
        
        m_manager->printPetscNotice(5,"Assemble the null space of PETSc matrix\n");
        m_manager->printPetscNotice(6,"Initialize the null space vectors\n");
        for (int i = 0; i < Dim; ++i)
        {
            MatCreateVecs(hessian,NULL,&nullvec[i]);
        }
        m_manager->printPetscNotice(6,"Set the values of the null space vectors\n");
        for (PetscInt i = Istart; i < Iend; ++i) 
        {
            //Instantiate the iterator at a particular position
            auto p_i = m_sysdata->particles.begin()+(int)(i/Dim);
            PetscInt id_i = abr::get<abr::id>(*p_i);
            setNullSpaceBasis(id_i, Istart, Iend, i);
        }
        m_manager->printPetscNotice(6,"Assemble and normalize \n");
        for (int i = 0; i < Dim; ++i)
        {
            VecAssemblyBegin(nullvec[i]);
            VecAssemblyEnd(nullvec[i]);
            VecNormalize(nullvec[i],NULL);
        }

        MatAssemblyBegin(hessian,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(hessian,MAT_FINAL_ASSEMBLY);
        MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_FALSE,2,nullvec,&constant);
        MatSetNullSpace(hessian,constant);
};

template< int Dim >
void PETScHessian<Dim>::setHessianValues(PetscInt id_i, PetscInt id_j, PetscInt Istart, PetscInt Iend, PetscInt real_id,Eigen::Matrix3d offdiag_ij)
{
    //Each "set value" must be carefully considered because we don't know the 
    //parallel layout of the matrix beforehand !! But we are allowerd to fill them row-by-row
    //row must also consider indexing at particular 
    //x-component of the row 
    if (Dim*id_i == real_id)
    { 
        //The main diagonal term
        MatSetValue(hessian,Dim*id_i, Dim*id_j, offdiag_ij(0,0),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i, Dim*id_j+1, offdiag_ij(0,1),ADD_VALUES);
        //The off diagonal term
        MatSetValue(hessian,Dim*id_i, Dim*id_i, -offdiag_ij(0,0),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i, Dim*id_i+1, -offdiag_ij(0,1),ADD_VALUES);
    }
    //y-component of the row
    else if (Dim*id_i+1 == real_id)
    {
        //the main diagonal term    
        MatSetValue(hessian,Dim*id_i+1, Dim*id_j, offdiag_ij(1,0),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i+1, Dim*id_j+1, offdiag_ij(1,1),ADD_VALUES);
        //The off diagonal term
        MatSetValue(hessian,Dim*id_i+1, Dim*id_i, -offdiag_ij(1,0),ADD_VALUES);
        MatSetValue(hessian,Dim*id_i+1, Dim*id_i+1, -offdiag_ij(1,1),ADD_VALUES);
    }
};

template< int Dim >
void PETScHessian<Dim>::setNullSpaceBasis(PetscInt id_i, PetscInt Istart, PetscInt Iend, PetscInt real_id)
{
    if (Dim*id_i == real_id)
    { 
        
        VecSetValue(nullvec[0],Dim*id_i,1, INSERT_VALUES);
        VecSetValue(nullvec[1],Dim*id_i,0, INSERT_VALUES);
        if (Dim == 3)
            VecSetValue(nullvec[Dim],Dim*id_i,0, INSERT_VALUES);
    }
    //y-component of the row
    else if (Dim*id_i+1 == real_id)
    {
        VecSetValue(nullvec[1],Dim*id_i+1,0, INSERT_VALUES);
        VecSetValue(nullvec[1],Dim*id_i+1,1, INSERT_VALUES);
        if (Dim == 3)
            VecSetValue(nullvec[Dim],Dim*id_i,0, INSERT_VALUES);
    }
};
/*
template< int Dim >
void PETScHessian<Dim>::setHessianValues(PetscInt id_i, PetscInt id_j, PetscInt Istart, PetscInt Iend, Eigen::Matrix3d offdiag_ij)
{
    //Each "set value" must be carefully considered because we don't know the 
    //parallel layout of the matrix beforehand !! But we are allowerd to fill them row-by-row
    //row must also consider indexing at particular 
    
    if (Dim*id_i >= Istart && Dim*id_i < Iend)
    {
        if (Dim*id_j >= Istart && Dim*id_j < Iend)
        {
            MatSetValue(hessian,Dim*id_i, Dim*id_j, offdiag_ij(0,0),ADD_VALUES);
        }
        if (Dim*id_j+1 >= Istart && Dim*id_j+1 < Iend)
        {
            MatSetValue(hessian,Dim*id_i, Dim*id_j+1, offdiag_ij(0,1),ADD_VALUES);
        }
        if (Dim == 3 && Dim*id_j+2 >= Istart && Dim*id_j+2 < Iend)
        {
            MatSetValue(hessian,Dim*id_i, Dim*id_j+2, offdiag_ij(0,2),ADD_VALUES);
        }
    }
    if (Dim*id_i+1 >= Istart && Dim*id_i+1 < Iend)
    {
        if (Dim*id_j >= Istart && Dim*id_j < Iend)
        {
            MatSetValue(hessian,Dim*id_i+1, Dim*id_j, offdiag_ij(1,0),ADD_VALUES);
        }
        if (Dim*id_j+1 >= Istart && Dim*id_j+1 < Iend)
        {
            MatSetValue(hessian,Dim*id_i+1, Dim*id_j+1, offdiag_ij(1,1),ADD_VALUES);
        }
        if (Dim == 3 && Dim*id_j+2 >= Istart && Dim*id_j+2 < Iend)
        {
            MatSetValue(hessian,Dim*id_i+1, Dim*id_j+2, offdiag_ij(1,2),ADD_VALUES);
        }
    }
    if (Dim == 3 && Dim*id_i+2 >= Istart && Dim*id_i+2 < Iend)
    {
        if (Dim*id_j >= Istart && Dim*id_j < Iend)
        {
            MatSetValue(hessian,Dim*id_i+2, Dim*id_j, offdiag_ij(2,0),ADD_VALUES);
        }
        if (Dim*id_j+1 >= Istart && Dim*id_j+1 < Iend)
        {
            MatSetValue(hessian,Dim*id_i+2, Dim*id_j+1, offdiag_ij(2,1),ADD_VALUES);
        }
        if (Dim == 3 && Dim*id_j+2 >= Istart && Dim*id_j+2 < Iend)
        {
            MatSetValue(hessian,Dim*id_i+2, Dim*id_j+2, offdiag_ij(2,2),ADD_VALUES);
        }
    }

    //The diagonal term 
    if (Dim*id_i >= Istart && Dim*id_i < Iend)
    {
        if (Dim*id_i >= Istart && Dim*id_i < Iend)
        {
            MatSetValue(hessian,Dim*id_i, Dim*id_i, -offdiag_ij(0,0),ADD_VALUES);
        }
        if (Dim*id_i+1 >= Istart && Dim*id_i+1 < Iend)
        {
            MatSetValue(hessian,Dim*id_i, Dim*id_i+1, -offdiag_ij(0,1),ADD_VALUES);
        }
        if (Dim == 3 && Dim*id_i+2 >= Istart && Dim*id_i+2 < Iend)
        {
            MatSetValue(hessian,Dim*id_i, Dim*id_i+2, -offdiag_ij(0,2),ADD_VALUES);
        }
    }
    if (Dim*id_i+1 >= Istart && Dim*id_i+1 < Iend)
    {
        if (Dim*id_i >= Istart && Dim*id_i < Iend)
        {
            MatSetValue(hessian,Dim*id_i+1, Dim*id_i, -offdiag_ij(1,0),ADD_VALUES);
        }
        if (Dim*id_i+1 >= Istart && Dim*id_i+1 < Iend)
        {
            MatSetValue(hessian,Dim*id_i+1, Dim*id_i+1, -offdiag_ij(1,1),ADD_VALUES);
        }
        if (Dim == 3 && Dim*id_i+2 >= Istart && Dim*id_i+2 < Iend)
        {
            MatSetValue(hessian,Dim*id_i+1, Dim*id_i+2, -offdiag_ij(1,2),ADD_VALUES);
        }
    }
    if (Dim == 3 && Dim*id_i+2 >= Istart && Dim*id_i+2 < Iend)
    {
        if (Dim*id_i >= Istart && Dim*id_i < Iend)
        {
            MatSetValue(hessian,Dim*id_i+2, Dim*id_i, -offdiag_ij(2,0),ADD_VALUES);
        }
        if (Dim*id_i+1 >= Istart && Dim*id_i+1 < Iend)
        {
            MatSetValue(hessian,Dim*id_i+2, Dim*id_i+1, -offdiag_ij(2,1),ADD_VALUES);
        }
        if (Dim == 3 && Dim*id_i+2 >= Istart && Dim*id_i+2 < Iend)
        {
            MatSetValue(hessian,Dim*id_i+2, Dim*id_i+2, -offdiag_ij(2,2),ADD_VALUES);
        }
    }
}
*/

template< class T >
void export_PETScHessian(py::module& m, const std::string& name)
{
    py::class_<T, PETScHessianBase, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr< PETScManager > , std::shared_ptr< MPI::Communicator > >())
    .def("destroyPETScObjects", &T::destroyPETScObjects)
    .def("assemblePETScObjects", &T::assemblePETScObjects)
    .def("setSystemData", &T::setSystemData)
    ;
};

#endif
