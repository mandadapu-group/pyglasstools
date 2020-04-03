#ifndef __BORN_TENSOR_H__
#define __BORN_TENSOR_H__

#include <pyglasstools/MathAndTypes.h>
#include <pyglasstools/potential/PairPotential.h>
namespace abr = Aboria;

class BornStiffnessTensor
{
    public:
        BornStiffnessTensor(int _dim) : dim(_dim){};
        ~BornStiffnessTensor(){};
        virtual Eigen::MatrixXd compute(const AboriaParticles::value_type& particle_i)
            {
                throw std::runtime_error("[ERROR] Born Stiffness Tensor needs information on the j-th particle and the potential");
                return Eigen::MatrixXd::Zero(dim,dim);
            }
        virtual Eigen::MatrixXd compute(const AboriaParticles::value_type& particle_i, 
                                        const AboriaParticles::value_type& particle_j)
            {
                throw std::runtime_error("[ERROR] Born Stiffness Tensor needs information on the potential force");
                return Eigen::MatrixXd::Zero(dim,dim);
            }
        
        virtual Eigen::MatrixXd compute(const AboriaParticles::value_type& particle_i, 
                                        const AboriaParticles::value_type& particle_j, 
                                        const std::shared_ptr<PairPotential>& potential)
            {
                
                //Compute pair-force!
                double forceval = potential->getPairForce();
                double stiffnessval = potential->getBondStiffness();
                Eigen::Vector3d F =  (forceval)*(potential->rij);
                Eigen::MatrixXd val =  Eigen::MatrixXd::Zero((int)3*(dim-1),(int)3*(dim-1));
                Eigen::Vector3d nij = potential->rij.normalized();

                //Compute virial stress
                if (dim > 2)
                {
                    //First row determines T_00 stress response
                    val(0,0) += (stiffnessval*(potential->rij[0]) + F[0])*potential->rij[0]*nij[0]*nij[0];
                    val(0,1) += (stiffnessval*(potential->rij[0]) + F[0])*potential->rij[0]*nij[1]*nij[1];
                    val(0,2) += (stiffnessval*(potential->rij[0]) + F[0])*potential->rij[0]*nij[0]*nij[1];
                    
                    //Second row determines T_11 stress response
                    val(1,0) += (stiffnessval*(potential->rij[1]) + F[1])*potential->rij[1]*nij[0]*nij[0];
                    val(1,1) += (stiffnessval*(potential->rij[1]) + F[1])*potential->rij[1]*nij[1]*nij[1];
                    val(1,2) += (stiffnessval*(potential->rij[1]) + F[1])*potential->rij[1]*nij[0]*nij[1];
                    
                    //Third row determines T_12 stress response
                    val(2,0) += (stiffnessval*(potential->rij[0]) + F[0])*potential->rij[1]*nij[0]*nij[0];
                    val(2,1) += (stiffnessval*(potential->rij[0]) + F[0])*potential->rij[1]*nij[1]*nij[1];
                    val(2,2) += (stiffnessval*(potential->rij[0]) + F[0])*potential->rij[1]*nij[0]*nij[1];
                    return val;
                }
                else
                {
                    throw std::runtime_error("No implementation in 3D yet");
                }
            }
        int dim;
};

#endif
