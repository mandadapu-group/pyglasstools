#ifndef __EIGEN_HESSIAN_BASE_H__
#define __EIGEN_HESSIAN_BASE_H__

#include <pyglasstools/elasticity/HessianBase.h>
#include "EigenManager.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/StdVector>

typedef Eigen::Triplet<double> Tripletd;
typedef Eigen::SparseMatrix<double> SparseMatd;

/*
 * Derived base class for constructing Hessian using Eigen's matrices and arrays.
 * It can be further specialized for different applications, e.g., normal mode analysis.
 */
class EigenHessianBase : public HessianBase
{
    public:
        
        std::shared_ptr< EigenManager > m_manager; //!< a Manager class for printing values from Eigen arrays and outputing error messages
        SparseMatd hessian; //!< the Hessian of the system, stored as a sparse Eigen matrix
        Eigen::MatrixXd misforce; //!< the mismatch force vector, denoted as Xi in many papers. 

        EigenHessianBase(   std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                            std::shared_ptr< EigenManager > manager);
        virtual ~EigenHessianBase()
        {
        }
        
};

/*
 * Parametrized constructor. Args:
 * sysdata: the particle system, storing all of relevant configurational data
 * potential: the pair potential of the chosen system.
 * manager: a manager class for printing Eigen arrays and error messages.
 *
 * All arguments are shared pointers, so that we are referring to an already-constructed objects. 
 */
EigenHessianBase::EigenHessianBase( std::shared_ptr< ParticleSystem > sysdata, std::shared_ptr< PairPotential > potential, 
                                    std::shared_ptr<EigenManager> manager)
    : HessianBase(sysdata,potential), m_manager(manager)
{
};

/*
 * Helper function to export EigenHessian base to Python
 */
void export_EigenHessianBase(py::module& m)
{
    py::class_<EigenHessianBase, HessianBase, std::shared_ptr<EigenHessianBase> >(m,"EigenHessianBase")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >, std::shared_ptr<EigenManager> >())
    ;
};

#endif
