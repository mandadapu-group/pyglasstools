#ifndef __PETSC_LINEAR_RESPONSE_H__
#define __PETSC_LINEAR_RESPONSE_H__

#include "PETScHessianBase.h"

class PYBIND11_EXPORT PETScLinearResponse
{
    protected:
        std::shared_ptr< PETScHessianBase > m_hessian; 
        
        //Properties of the particle pair
        double forcedipole;
        PetscInt id_i, id_j;
        Eigen::Vector3d rij; 
        
        //List of observables to be computed!
        std::map<std::string, std::shared_ptr< PETScVectorFieldBase > > m_vectorfields;
        std::map<std::string, std::shared_ptr< GlobalPropertyBase > > m_observables;
        
        //PETSc error code 
        PetscErrorCode ierr; 
    public:
        PETScLinearResponse( std::shared_ptr< PETScHessianBase > hessian)
            : m_hessian(hessian), forcedipole(std::numeric_limits<double>::max()), id_i(0), id_j(0)
        {
            //Add any relevant command-line options for the PETSc linear algebra solver.
            ierr = PetscOptionsInsertString(NULL,m_hessian->m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
        }
        virtual ~PETScLinearResponse()
        {
        };
        
        void setHessian(std::shared_ptr< PETScHessianBase > hessian)
        {
            m_hessian = hessian;
        }
        
        /* Add a new vector field observable to a list of existing ones. Args:
         * obs: the new vector field observable
         */
        virtual void addVectorField(const std::shared_ptr< PETScVectorFieldBase >& obs)
        {
            m_vectorfields.insert(std::pair<std::string, std::shared_ptr< PETScVectorFieldBase > >(obs->name, obs));
        }
        
        /* Add a new global observable/property toa  list of existing ones. Args:
         * obs: the new global observable
         */
        virtual void addGlobalProperty(const std::shared_ptr< GlobalPropertyBase >& obs)
        {
            m_observables.insert(std::pair<std::string, std::shared_ptr< GlobalPropertyBase > >(obs->name, obs));
        }
        
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
        
        void solveLinearResponseProblem();
        void findRandomPairForce();
        void findMinimumPairForce();
        void setLocLandscapeForce(const Vec& forcing, std::string mode);
        //void saveForceDipoleProblem(std::shared_ptr<MPI::LogFile> logfile);
};

void PETScLinearResponse::solveLinearResponseProblem()
{
    if (m_vectorfields.count("displacement") > 0)
    {
        KSP ksp;
        Vec forcing;
        
        MatCreateVecs(m_hessian->hessian,NULL,&forcing);
        
        //Construct the KSP object with hessian as our matrix
        KSPCreate(PETSC_COMM_WORLD,&ksp);
        KSPSetOperators(ksp,m_hessian->hessian,m_hessian->hessian);
        
        setForceVector(forcing);
        KSPSetFromOptions(ksp);
        
        //Do another insert command line options just to be sure
        //m_hessian->ierr = PetscOptionsInsertString(NULL,m_hessian->m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,m_hessian->ierr);
        //Create the vector to store solution
        m_vectorfields["displacement"]->createVector(m_hessian->hessian);
        
        //Save pair particle properties chosen by 
        Eigen::Vector3d center;
        center <<  0.5*(2*abr::get<position>(m_hessian->m_sysdata->particles[id_i])[0]-rij[0]),
                   0.5*(2*abr::get<position>(m_hessian->m_sysdata->particles[id_i])[1]-rij[1]),
                   0.5*(2*abr::get<position>(m_hessian->m_sysdata->particles[id_i])[2]-rij[2]);
        center = m_hessian->m_sysdata->simbox->minImage(center);
        for (int i = 0; i < center.size(); ++i)
        {
            if (m_observables.count("/forcedipole_center") > 0) m_observables["forcedipole_center"]->setValue(center[i],i);
            if (m_observables.count("forcedipole_pi") > 0) m_observables["forcedipole_pi"]->setValue(abr::get<position>(m_hessian->m_sysdata->particles[id_i])[i],i);
            if (m_observables.count("forcedipole_pj") > 0) m_observables["forcedipole_pj"]->setValue(abr::get<position>(m_hessian->m_sysdata->particles[id_i])[i]-rij[i],i);
        }
        
        m_hessian->m_manager->printPetscNotice(5, "Solve the Hu=f problem. \n");
        KSPSolve(ksp,forcing,m_vectorfields["displacement"]->vectorobs);
        m_hessian->m_manager->printPetscNotice(5, "Calculation finished. \n");
        
        KSPDestroy(&ksp);
        VecDestroy(&forcing);
    }
};

void PETScLinearResponse::findRandomPairForce()
{
    forcedipole = std::numeric_limits<double>::max();
    if (m_hessian->m_comm->isRoot())
    {
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist6(0,m_hessian->m_sysdata->particles.size()-1);
        std::stringstream string_stream; string_stream << "Looking for pairs of particles to perturb" << std::endl; 
        m_hessian->m_manager->printPetscNotice(5,string_stream.str());
        //for( auto p_i = m_hessian->m_sysdata->particles.begin(); p_i != m_hessian->m_sysdata->particles.end(); ++p_i)
        bool notfound = true;
        do
        {
            //std::cout << dist6(rng) << std::endl;
            //py::print(id_i,id_i+1,Iend);
            auto p_i = m_hessian->m_sysdata->particles.begin()+dist6(rng);
            PetscInt tempid_i = abr::get<abr::id>(*p_i);
            for( auto p_j = abr::euclidean_search(m_hessian->m_sysdata->particles.get_query(), 
                        abr::get<position>(*p_i), m_hessian->max_rcut); p_j != false; ++p_j)
            {
                PetscInt tempid_j = abr::get<abr::id>(*p_j);
                 
                //This needs to be changed so that it works in 2/3 dimensions 
                if (tempid_i != tempid_j)
                {
                    //Compute a list of local obsercavles 
                    //Next we loop through the j-th particles for the virial stress
                    //set the distance between particle i and particle j
                    Eigen::Vector3d bar_rij(-p_j.dx()[0], -p_j.dx()[1], -p_j.dx()[2]); //this is rij = ri-rj, with PBC taken into account
                    
                    //Don't forget to set diameters of the potential
                    double di =  abr::get<diameter>(*p_i);
                    double dj =  abr::get<diameter>(*p_j);
                    
                    //Make sure the particle is unique
                    //if (m_hessian->m_potential->getRcut() > rij.dot(rij) && abs(m_hessian->m_potential->getPairForceDivR()) < abs(force_threshold))
                    double testdipole = abs(m_hessian->m_potential->getPairForceDivR(bar_rij,di,dj));
                    if (    m_hessian->m_potential->getRcut(bar_rij,di,dj) > bar_rij.dot(bar_rij) && 
                            testdipole*bar_rij.norm() > m_hessian->m_manager->fd_random_min 
                            && testdipole*bar_rij.norm() < m_hessian->m_manager->fd_random_max)
                    {
                        string_stream.str(std::string());
                        string_stream.clear();
                        string_stream << "Found it!" << std::endl;
                        m_hessian->m_manager->printPetscNotice(5,string_stream.str()); 
                        
                        string_stream.str(std::string());
                        string_stream.clear();
                        string_stream << "Force min: " << m_hessian->m_manager->fd_random_min << " and Force max: " << m_hessian->m_manager->fd_random_max << " distance is " << bar_rij[0] << " " << bar_rij[1] << std::endl;
                        m_hessian->m_manager->printPetscNotice(5,string_stream.str()); 
                        
                        forcedipole = 1.0/bar_rij.norm();//m_hessian->m_potential->getPairForceDivR(bar_rij,di,dj);
                        rij = bar_rij;
                        id_i = tempid_i; 
                        id_j = tempid_j;
                        notfound = false; 
                        break;
                    }
                }
            }
        }
        while (notfound);
    }
    //Now, broadcast the results for id_i, id_j, rij, and forcedipole
    //Now determine which rank has the smallest force. Once determined, we will broadcast the values
    m_hessian->m_comm->bcast(id_i,0); 
    m_hessian->m_comm->bcast(id_j,0); 
    m_hessian->m_comm->bcast(forcedipole,0); 
    m_hessian->m_comm->bcast(rij,0);
}

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

void PETScLinearResponse::setLocLandscapeForce(const Vec& forcing, std::string mode)
{
    //All processes will look for
    //Now determine which rank has the smallest force. Once determined, we will broadcast the values
    forcedipole = 0.0;//std::numeric_limits<double>::max();
    m_hessian->m_manager->notice(5) << "Looking for pairs of particles to perturb" << std::endl; 
    PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(m_hessian->hessian, &Istart, &Iend);
    //for( auto p_i = m_hessian->m_sysdata->particles.begin(); p_i != m_hessian->m_sysdata->particles.end(); ++p_i)
    for (PetscInt i = Istart; i < Iend; ++i) 
    {
        //Instantiate the iterator at a particular position
        auto p_i = m_hessian->m_sysdata->particles.begin()+(int)(i/2); //Divide by i/Dim to get the actual particle index
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
                if (m_hessian->m_potential->getRcut(bar_rij,di,dj) > bar_rij.dot(bar_rij))// && abs(forcedipole) < abs(testdipole) && testdipole > 0 )
                {
                    double factor = 0;
                    if (mode == "nonaffine")
                    {
                        factor = m_hessian->m_potential->getBondStiffness(bar_rij, di, dj)+m_hessian->m_potential->getPairForceDivR(bar_rij, di, dj);
                        factor *= bar_rij[0]*bar_rij[1]/(bar_rij.norm()*bar_rij.norm());
                        //double factor1 = 0.25*(bar_rij[0]*bar_rij[0]-bar_rij[1]*bar_rij[1])*(bar_rij[0]*bar_rij[0]-bar_rij[1]*bar_rij[1]);
                        //double factor2 = bar_rij[0]*bar_rij[0]*bar_rij[1]*bar_rij[1];
                        //factor *= -sqrt(factor1+factor2)/(bar_rij.norm()*bar_rij.norm());
                    }
                    else if(mode == "uniform")
                    {
                        factor = 1.0;
                    }
                    
                    if (2*tempid_i == i)
                    { 
                        VecSetValue(forcing,2*tempid_i,factor*bar_rij[0], ADD_VALUES);
                    }
                    else if (2*tempid_i+1 ==  i)
                    {
                        VecSetValue(forcing,2*tempid_i+1,factor*bar_rij[1], ADD_VALUES);
                    }
                }
            }
        }
    }
}

void export_PETScLinearResponse(py::module& m)
{
    py::class_< PETScLinearResponse, std::shared_ptr< PETScLinearResponse > >(m,"PETScLinearResponse")
    .def(py::init< std::shared_ptr< PETScHessianBase > >())
    .def("setHessian", &PETScLinearResponse::setHessian)
    .def("solveLinearResponseProblem", &PETScLinearResponse::solveLinearResponseProblem)
    .def("addVectorField", &PETScLinearResponse::addVectorField)
    .def("addGlobalProperty", &PETScLinearResponse::addGlobalProperty)
    ;
};


#endif //__PETSC_LINEAR_RESPONSE_H__
