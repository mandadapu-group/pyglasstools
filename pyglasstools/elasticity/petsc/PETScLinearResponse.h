#ifndef __PETSC_LINEAR_RESPONSE_H__
#define __PETSC_LINEAR_RESPONSE_H__

#include "PETScHessianBase.h"

class PYBIND11_EXPORT PETScLinearResponse : public PETScCalculatorBase
{
    protected:
        //Properties of the particle pair
        //double forcedipole;
        //PetscInt id_i, id_j;
        //Eigen::Vector3d rij; 
        
        //The list of \Xi vectors
        unsigned int xivecsdim; 
        Vec* xivectors;

        //The KSP Objects
        KSP ksp;
    public:
        PETScLinearResponse( std::shared_ptr< PETScHessianBase > hessian)
            : PETScCalculatorBase(hessian)
        {
            int Dim = m_hessian->m_sysdata->simbox->dim;
            xivecsdim = (unsigned int)Dim*((unsigned int)Dim+1)/2;
            assembleObjects();     
        }
        virtual ~PETScLinearResponse()
        {
            destroyObjects();
        }
        
        /* Build the \Xi_{ijk}^\gamma vector
         * needed for computing non-affine elasticity tensor
         */ 
        void buildXiVector();
        
        void destroyObjects()
        {
            for(unsigned int i = 0; i < xivecsdim; ++i)
            {
                VecDestroy(&xivectors[i]);
            }
            PetscFree(xivectors);
            KSPDestroy(&ksp);
        } 
        void assembleObjects()
        {
            buildXiVector();
            //Construct the KSP object with hessian as our matrix
            KSPCreate(PETSC_COMM_WORLD,&ksp);
            KSPSetOperators(ksp,m_hessian->hessian,m_hessian->hessian);
            KSPSetFromOptions(ksp);
        }        
        /* Computing the non-affine part of the elasticity tensor. Args:
         * eps: class holding all eigenmodes.
         */
        void calculateNonAffineTensor();
        
        /* Finding a pair of particles to apply force dipoles on
         *
         */ 
        //void findMinimumPairForce();
        
        /*
         * Helper function to set the component of the mismatch force vector, denoted as Xi in many papers.
         */ 
        void setXiVectorValues(PetscInt id_i,PetscInt real_id,double factor,Eigen::Vector3d rij, Eigen::Vector3d nij);
        
        //void saveForceDipoleProblem(std::shared_ptr<MPI::LogFile> logfile);
        //virtual void setForceVector(const Vec& forcing)
};

/* 
virtual void setForceVector(const Vec& forcing)
{
    std::stringstream streamm;
    streamm << "Assembling forcing vector" << std::endl;
    m_hessian->m_manager->printPetscNotice(5,streamm.str());
     
    PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(m_hessian->hessian, &Istart, &Iend);
    if (m_hessian->m_manager->fd_mode == "random")
    {
        findRandomPairForce();
        
        streamm << "Perturbing particle " << id_i << " and particle " << id_j << " with force_x " << forcedipole*rij[0] << " with force_y " << forcedipole*rij[1] << std::endl; 
        m_hessian->m_manager->printPetscNotice(5, streamm.str());
        streamm.str(std::string());
        streamm.clear();

        if (2*id_i < Iend && 2*id_i >= Istart)
        {
            VecSetValue(forcing,2*id_i,-forcedipole*rij[0], ADD_VALUES);
            VecSetValue(forcing,2*id_i+1,-forcedipole*rij[1], ADD_VALUES);
        }
        if (2*id_j < Iend && 2*id_j >= Istart)
        {
            VecSetValue(forcing,2*id_j,forcedipole*rij[0], ADD_VALUES);
            VecSetValue(forcing,2*id_j+1,forcedipole*rij[1], ADD_VALUES);
        }
    }
    else if (m_hessian->m_manager->fd_mode == "uniform")
    {
        setLocLandscapeForce(forcing,"uniform");
        streamm << "Applying force dipoles with uniform magnitude to every possible pairs" << std::endl;
        m_hessian->m_manager->printPetscNotice(5, streamm.str());
        streamm.str(std::string());
        streamm.clear();
    }
    else if (m_hessian->m_manager->fd_mode == "nonaffine")
    {
        setLocLandscapeForce(forcing,"nonaffine");
        streamm << "Applying forces that induce nonaffine displacement" << std::endl;
        m_hessian->m_manager->printPetscNotice(5, streamm.str());
        streamm.str(std::string());
        streamm.clear();
    }
    else
    {
        streamm << "Option not available! Only 'random', 'uniform', or 'nonaffine' are available" << std::endl;
        m_hessian->m_manager->printPetscNotice(5, streamm.str());
        streamm.str(std::string());
        streamm.clear();
    }
    std::stringstream string_stream;
    string_stream << "Assembling forcing vector" << std::endl; 
    m_hessian->m_manager->printPetscNotice(5,string_stream.str());
    VecAssemblyBegin(forcing);
    VecAssemblyEnd(forcing);
}
*/

    

/*
void PETScLinearResponse::findMinimumPairForce()
{
    //All processes will look for
    //Now determine which rank has the smallest force. Once determined, we will broadcast the values
    forcedipole = 0.0;//std::numeric_limits<double>::max();
    if(m_hessian->m_comm->isRoot())
    {
        m_hessian->m_manager->notice(5) << "Looking for pairs of particles to perturb" << std::endl; 
        for( auto p_i = m_hessian->m_sysdata->particles.begin(); p_i != m_hessian->m_sysdata->particles.end(); ++p_i)
        {
            //std::cout << dist6(rng) << std::endl;
            //py::print(id_i,id_i+1,Iend);
            for( auto p_j = abr::euclidean_search(m_hessian->m_sysdata->particles.get_query(), 
                        abr::get<position>(*p_i), m_hessian->max_rcut); p_j != false; ++p_j)
            {
                PetscInt tempid_j = abr::get<abr::id>(*p_j);
                PetscInt tempid_i = abr::get<abr::id>(*p_i);
               
                //This needs to be changed so that it works in 2/3 dimensions 
                if (tempid_i != tempid_j)
                {
                    //Compute a list of local obsercavles 
                    //Next we loop through the j-th particles for the virial stress
                    //set the distance between particle i and particle j
                    Eigen::Vector3d bar_rij(-p_j.dx()[0], -p_j.dx()[1], -p_j.dx()[2]);
                    
                    //Don't forget to set diameters of the potential
                    double di =  abr::get<diameter>(*p_i);
                    double dj =  abr::get<diameter>(*p_j);
                    
                    //Make sure the particle is unique
                    //if (m_hessian->m_potential->getRcut() > rij.dot(rij) && abs(m_hessian->m_potential->getPairForceDivR()) < abs(force_threshold))
                    double testdipole = m_hessian->m_potential->getPairForceDivR(bar_rij,di,dj)*bar_rij.norm();
                    if (m_hessian->m_potential->getRcut(bar_rij,di,dj) > bar_rij.dot(bar_rij) && abs(forcedipole) < abs(testdipole) && testdipole > 0 )
                    {
                        forcedipole = m_hessian->m_potential->getPairForceDivR(bar_rij,di,dj);
                        rij = bar_rij;
                        id_i = tempid_i; 
                        id_j = tempid_j; 
                    }
                }
            }
        }
    }
    m_hessian->m_comm->barrier();
    m_hessian->m_comm->bcast(id_i,0); 
    m_hessian->m_comm->bcast(id_j,0); 
    m_hessian->m_comm->bcast(forcedipole,0); 
    m_hessian->m_comm->bcast(rij,0);
}
*/

void PETScLinearResponse::buildXiVector()
{
    int Dim = m_hessian->m_sysdata->simbox->dim;
    //Compute the maximum cut-off radius. This is set by either the force cut-off radius or the diameter
    double maxforce_rcut = m_hessian->m_potential->scaled_rcut;
    double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_hessian->m_sysdata->particles)), 
                                            std::end(abr::get<diameter>(m_hessian->m_sysdata->particles)) );
    maxforce_rcut *= maxdiameter;
    double max_rcut = std::max(maxforce_rcut,maxdiameter); 
    
    //Construct the PETSc hessian matrix, but don't assemble yet!
    PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(m_hessian->hessian, &Istart, &Iend);
    m_hessian->m_manager->printPetscNotice(5,"Assembling the Xi_{ijkl}^gamma vector \n");
    PetscMalloc1(xivecsdim,&xivectors);    
    for (unsigned int i = 0; i < xivecsdim; ++i)
    {
        MatCreateVecs(m_hessian->hessian,NULL,&xivectors[i]);
    }
    
    //The first-loop is over the particle index and Cartesian component, as owned by the processor
    for (PetscInt i = Istart; i < Iend; ++i) 
    {
        //Instantiate the iterator at a particular position
        auto p_i = m_hessian->m_sysdata->particles.begin()+(int)(i/Dim); //Divide by i/Dim to get the actual particle index
        
        //The second for-loop is to go through all possible neighbors of the i-th particle
        for( auto p_j = abr::euclidean_search(m_hessian->m_sysdata->particles.get_query(), 
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
                if (m_hessian->m_potential->getRcut(rij, di, dj) > rij.dot(rij)) 
                {
                    //Note PairForceDivR is -phi_r(r)/r, hence the negative sign is already encoded inside factor
                    
                    Eigen::Matrix3d offdiag_ij;
                    double factor;
                    factor = m_hessian->m_potential->getBondStiffness(rij, di, dj)+m_hessian->m_potential->getPairForceDivR(rij, di, dj);
                    offdiag_ij = -factor*nij*nij.transpose()
                                                 +Eigen::Matrix3d::Identity()*m_hessian->m_potential->getPairForceDivR(rij, di, dj);
                    
                    setXiVectorValues(id_i,i,factor,rij,nij);
                }
            }
        }
    }
    
    for (unsigned int i = 0; i < xivecsdim; ++i)
    {
        VecAssemblyBegin(xivectors[i]);
        VecAssemblyEnd(xivectors[i]);
    }
};

void PETScLinearResponse::setXiVectorValues(PetscInt id_i, int real_id, double factor, Eigen::Vector3d rij, Eigen::Vector3d nij)
{
    int Dim = m_hessian->m_sysdata->simbox->dim;
    if (Dim*id_i == real_id)
    {

        VecSetValue(xivectors[0],Dim*id_i,-factor*rij[0]*nij[0]*nij[0],ADD_VALUES); 
        VecSetValue(xivectors[1],Dim*id_i,-factor*rij[0]*nij[1]*nij[0],ADD_VALUES); 
        VecSetValue(xivectors[2],Dim*id_i,-factor*rij[1]*nij[1]*nij[0],ADD_VALUES); 

        if (Dim == 3)
        {
            VecSetValue(xivectors[3],Dim*id_i,-factor*rij[1]*nij[2]*nij[0],ADD_VALUES); 
            VecSetValue(xivectors[4],Dim*id_i,-factor*rij[2]*nij[2]*nij[0],ADD_VALUES); 
            VecSetValue(xivectors[5],Dim*id_i,-factor*rij[0]*nij[2]*nij[0],ADD_VALUES); 
        }
    }
    
    //y-component of the row
    else if (Dim*id_i+1 == real_id)
    {
        VecSetValue(xivectors[0],Dim*id_i+1,-factor*rij[0]*nij[0]*nij[1],ADD_VALUES); 
        VecSetValue(xivectors[1],Dim*id_i+1,-factor*rij[0]*nij[1]*nij[1],ADD_VALUES); 
        VecSetValue(xivectors[2],Dim*id_i+1,-factor*rij[1]*nij[1]*nij[1],ADD_VALUES); 
        if (Dim == 3)
        {
            VecSetValue(xivectors[3],Dim*id_i+1,-factor*rij[1]*nij[2]*nij[1],ADD_VALUES); 
            VecSetValue(xivectors[4],Dim*id_i+1,-factor*rij[2]*nij[2]*nij[1],ADD_VALUES); 
            VecSetValue(xivectors[5],Dim*id_i+1,-factor*rij[0]*nij[2]*nij[1],ADD_VALUES); 
        }
    }

    //z-component of the row
    else if (Dim*id_i+2 == real_id && Dim == 3)
    {
        VecSetValue(xivectors[0],Dim*id_i+2,-factor*rij[0]*nij[0]*nij[2],ADD_VALUES); 
        VecSetValue(xivectors[1],Dim*id_i+2,-factor*rij[0]*nij[1]*nij[2],ADD_VALUES); 
        VecSetValue(xivectors[2],Dim*id_i+2,-factor*rij[1]*nij[1]*nij[2],ADD_VALUES); 
        VecSetValue(xivectors[3],Dim*id_i+2,-factor*rij[1]*nij[2]*nij[2],ADD_VALUES); 
        VecSetValue(xivectors[4],Dim*id_i+2,-factor*rij[2]*nij[2]*nij[2],ADD_VALUES); 
        VecSetValue(xivectors[5],Dim*id_i+2,-factor*rij[0]*nij[2]*nij[2],ADD_VALUES); 
    }
    
}

void PETScLinearResponse::calculateNonAffineTensor()//const EPS& eps) 
{
    int Dim = m_hessian->m_sysdata->simbox->dim;
    if (m_observables.count("nonaffinetensor") && m_hessian->areDiagonalsNonZero())
    {
        m_hessian->m_manager->printPetscNotice(5,"Computing the nonaffine elasticity tensor \n");
        
        m_observables["nonaffinetensor"]->clear();
        Vec displacement;
        
        //We need a bit of prep because the layout may or may not be compatible between hessian and misforce
        MatCreateVecs(m_hessian->hessian,NULL,&displacement);
        
        //Although the non-affine tensor is stores as a proper 4-th rank tensor, the mismatch force vector, i.e., misforce
        //is stored as a flattened array. There is a specific mapping between the flattened array to the 4-th rank tensor. 
        for(int l = 0; l < m_observables["nonaffinetensor"]->getDimension(3); ++l)
        {
            for(int k = 0; k < m_observables["nonaffinetensor"]->getDimension(2); ++k)
            {
                for(int j = 0; j < m_observables["nonaffinetensor"]->getDimension(1); ++j)
                {
                    for (int i = 0; i < m_observables["nonaffinetensor"]->getDimension(0); ++i)
                    {

                        //if i,j or k,l  is 0,0 ---> 0
                        //if i,j or k,l  is 0,1 ---> 1
                        //if i,j or k,l  is 1,1 ---> 1
                        PetscScalar val = 0;
                        int index = l+Dim*(k+Dim*(j+Dim*i));
                        int id_i = i+j;
                        int id_j = k+l;
                        if (Dim == 3 && (id_i == 2 || id_j == 2))
                        {
                            if ((i == 0 && j == 2) || (i == 2 && j == 0) )
                            {
                                id_i = 5;
                            }
                            if ((k == 0 && l == 2) || (k == 2 && l == 0))
                            {
                                id_j = 5;
                            }
                            KSPSolve(ksp,xivectors[id_j],displacement);
                            VecDot(xivectors[id_i],displacement,&val);
                            m_observables["nonaffinetensor"]->addValue(val,index);
                        }
                        else
                        {
                            KSPSolve(ksp,xivectors[id_j],displacement);
                            VecDot(xivectors[id_i],displacement,&val);
                            m_observables["nonaffinetensor"]->addValue(val,index);
                        }
                    }
                }
            }
        }
        m_observables["nonaffinetensor"]->multiplyBy(1/m_hessian->m_sysdata->simbox->vol);
        VecDestroy(&displacement);
    };
};

void export_PETScLinearResponse(py::module& m)
{
    py::class_< PETScLinearResponse, std::shared_ptr< PETScLinearResponse > >(m,"PETScLinearResponse")
    .def(py::init< std::shared_ptr< PETScHessianBase > >())
    .def("calculateNonAffineTensor", &PETScLinearResponse::calculateNonAffineTensor)
    .def("setHessian", &PETScLinearResponse::setHessian)
    .def("addVectorField", &PETScLinearResponse::addVectorField)
    .def("addGlobalProperty", &PETScLinearResponse::addGlobalProperty)
    .def("assembleObjects", &PETScLinearResponse::assembleObjects)
    .def("destroyObjects", &PETScLinearResponse::destroyObjects)
    ;
};


#endif //__PETSC_LINEAR_RESPONSE_H__
