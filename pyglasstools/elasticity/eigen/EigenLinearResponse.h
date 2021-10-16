#ifndef __EIGEN_LINEAR_RESPONSE_H__
#define __EIGEN_LINEAR_RESPONSE_H__

#include "EigenHessianBase.h"
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
#include "LSQR.h"

class PYBIND11_EXPORT EigenLinearResponse
{
    protected:
        std::shared_ptr< EigenHessianBase > m_hessian; 
        
        //Properties of the particle pair
        double forcedipole;
        int id_i, id_j;
        Eigen::Vector3d rij; 
        
        //List of observables to be computed!
        std::map<std::string, std::shared_ptr< EigenVectorFieldBase > > m_vectorfields;
        std::map<std::string, std::shared_ptr< GlobalPropertyBase > > m_observables;
        
    public:
        EigenLinearResponse( std::shared_ptr< EigenHessianBase > hessian)
            : m_hessian(hessian), forcedipole(std::numeric_limits<double>::max()), id_i(0), id_j(0)
        {
        }
        virtual ~EigenLinearResponse()
        {
        };
        
        void setHessian(std::shared_ptr< EigenHessianBase > hessian)
        {
            m_hessian = hessian;
        }
        
        /* Add a new vector field observable to a list of existing ones. Args:
         * obs: the new vector field observable
         */
        virtual void addVectorField(const std::shared_ptr< EigenVectorFieldBase >& obs)
        {
            m_vectorfields.insert(std::pair<std::string, std::shared_ptr< EigenVectorFieldBase > >(obs->name, obs));
        }
        
        /* Add a new global observable/property toa  list of existing ones. Args:
         * obs: the new global observable
         */
        virtual void addGlobalProperty(const std::shared_ptr< GlobalPropertyBase >& obs)
        {
            m_observables.insert(std::pair<std::string, std::shared_ptr< GlobalPropertyBase > >(obs->name, obs));
        }
        
        void solveLinearResponseProblem();
        void setLocLandscapeForce(Eigen::VectorXd& forcing);
        void findRandomPairForce();
        void findMinimumPairForce();
        
        //void saveForceDipoleProblem(std::shared_ptr<MPI::LogFile> logfile);
        virtual void setForceVector(Eigen::VectorXd& forcing)
        {
            
            findMinimumPairForce();
            //findRandomPairForce();
            m_hessian->m_manager->notice(5) << "Assembling forcing vector" << std::endl;
            std::stringstream streamm;
            streamm << "Perturbing particle " << id_i << " and particle " << id_j << " with force_x " << forcedipole*rij[0] << " with force_y " << forcedipole*rij[1] << std::endl; 
            m_hessian->m_manager->notice(5) << streamm.str();
            forcing[2*id_i] = -forcedipole*rij[0];
            forcing[2*id_i+1] = -forcedipole*rij[1];
            forcing[2*id_j] = forcedipole*rij[0];
            forcing[2*id_j+1] = forcedipole*rij[1];
            
        }
};

void EigenLinearResponse::solveLinearResponseProblem()
{
    if (m_vectorfields.count("displacement") > 0)
    {
        Eigen::VectorXd forcing(m_hessian->hessian.rows());
        setLocLandscapeForce(forcing);
        //setForceVector(forcing);
        //Create the vector to store solution
        //Eigen::SPQR < SparseMatrix<double> > solver;
        //Eigen::LeastSquaresConjugateGradient < SparseMatrix<double> > solver;
        //Eigen::BiCGSTAB < SparseMatrix<double> > solver;
       	//Eigen::SparseLU < SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
        Eigen::SparseQR < SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
        m_hessian->hessian.makeCompressed();
        solver.compute(m_hessian->hessian);
        // decomposition failed
        if (solver.info() != Eigen::Success)
        {
              std::cout << "Decomposition Failed" << std::endl; 
        }
        m_hessian->m_manager->notice(5) << "Solve the Hu=f problem. " << std::endl;
        Eigen::VectorXd xcg = solver.solve(forcing);
        if (solver.info() != Eigen::Success)
        {
                std::cout << "Solving Failed" << std::endl;
        }
		/*
        */
		/*
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver(m_hessian->hessian);
        if (solver.info() != Eigen::Success)
                std::cout << "CG::Decomposition Failed" << std::endl;
        Eigen::VectorXd xcg = solver.solve(forcing);
        if (solver.info() != Eigen::Success)
                std::cout << "CG::Solving Failed" << std::endl;
		*/
        m_vectorfields["displacement"]->vectorobs = xcg;
        m_hessian->m_manager->notice(5) << "Calculation finished." << std::endl;
        
    }
};

void EigenLinearResponse::findRandomPairForce()
{
    forcedipole = std::numeric_limits<double>::max();
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(0,m_hessian->m_sysdata->particles.size()-1);
    m_hessian->m_manager->notice(5) << "Looking for pairs of particles to perturb" << std::endl; 
    bool notfound = true;
    do
    {
        auto p_i = m_hessian->m_sysdata->particles.begin()+dist6(rng);
        int tempid_i = abr::get<abr::id>(*p_i);
        for( auto p_j = abr::euclidean_search(m_hessian->m_sysdata->particles.get_query(), 
                    abr::get<position>(*p_i), m_hessian->max_rcut); p_j != false; ++p_j)
        {
            int tempid_j = abr::get<abr::id>(*p_j);
             
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
                    m_hessian->m_manager->notice(5) << "Found it!" << std::endl;
                    m_hessian->m_manager->notice(5) << "Force min: " << m_hessian->m_manager->fd_random_min << " and Force max: " << m_hessian->m_manager->fd_random_max << " distance is " << bar_rij[0] << " " << bar_rij[1] << std::endl;
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

void EigenLinearResponse::findMinimumPairForce()
{
    //All processes will look for
    //Now determine which rank has the smallest force. Once determined, we will broadcast the values
    forcedipole = 0.0;//std::numeric_limits<double>::max();
    m_hessian->m_manager->notice(5) << "Looking for pairs of particles to perturb" << std::endl; 
    for( auto p_i = m_hessian->m_sysdata->particles.begin(); p_i != m_hessian->m_sysdata->particles.end(); ++p_i)
    {
        //std::cout << dist6(rng) << std::endl;
        //py::print(id_i,id_i+1,Iend);
        for( auto p_j = abr::euclidean_search(m_hessian->m_sysdata->particles.get_query(), 
                    abr::get<position>(*p_i), m_hessian->max_rcut); p_j != false; ++p_j)
        {
            int tempid_j = abr::get<abr::id>(*p_j);
            int tempid_i = abr::get<abr::id>(*p_i);
           
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

void EigenLinearResponse::setLocLandscapeForce(Eigen::VectorXd& forcing)
{
    //All processes will look for
    //Now determine which rank has the smallest force. Once determined, we will broadcast the values
    forcedipole = 0.0;//std::numeric_limits<double>::max();
    m_hessian->m_manager->notice(5) << "Looking for pairs of particles to perturb" << std::endl; 
    unsigned int hessian_length = (unsigned int)2*m_hessian->m_sysdata->particles.size();
    //for( auto p_i = m_hessian->m_sysdata->particles.begin(); p_i != m_hessian->m_sysdata->particles.end(); ++p_i)
    double totforcex = 0;
    double totforcey = 0;
    for (unsigned int i = 0; i < hessian_length; ++i) 
    {
        //Instantiate the iterator at a particular position
        auto p_i = m_hessian->m_sysdata->particles.begin()+(int)(i/2); //Divide by i/Dim to get the actual particle index
        //std::cout << dist6(rng) << std::endl;
        //py::print(id_i,id_i+1,Iend);
        for( auto p_j = abr::euclidean_search(m_hessian->m_sysdata->particles.get_query(), 
                    abr::get<position>(*p_i), m_hessian->max_rcut); p_j != false; ++p_j)
        {
            int tempid_j = abr::get<abr::id>(*p_j);
            int tempid_i = abr::get<abr::id>(*p_i);
           
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
                    //double testdipole = m_hessian->m_potential->getPairForceDivR(bar_rij,di,dj)*bar_rij.norm();
                    double factor = m_hessian->m_potential->getBondStiffness(bar_rij, di, dj)+m_hessian->m_potential->getPairForceDivR(bar_rij, di, dj);
                    if (2*tempid_i == i)
                    { 
                        forcing[2*tempid_i] += factor*bar_rij[0];//, ADD_VALUES);
                        totforcex += factor*bar_rij[0];
                    }
                    else if (2*tempid_i+1 ==  i)
                    {
                        forcing[2*tempid_i+1] += factor*bar_rij[1];//, ADD_VALUES);
                        totforcey += factor*bar_rij[1];
                    }//forcedipole = m_hessian->m_potential->getPairForceDivR(bar_rij,di,dj);
                    //rij = bar_rij;
                    //id_i = tempid_i; 
                    //id_j = tempid_j; 
                }
            }
        }
    }
}

void export_EigenLinearResponse(py::module& m)
{
    py::class_< EigenLinearResponse, std::shared_ptr< EigenLinearResponse > >(m,"EigenLinearResponse")
    .def(py::init< std::shared_ptr< EigenHessianBase > >())
    .def("setHessian", &EigenLinearResponse::setHessian)
    .def("solveLinearResponseProblem", &EigenLinearResponse::solveLinearResponseProblem)
    .def("addVectorField", &EigenLinearResponse::addVectorField)
    .def("addGlobalProperty", &EigenLinearResponse::addGlobalProperty)
    ;
};


#endif //__EIGEN_LINEAR_RESPONSE_H__
