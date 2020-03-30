#ifndef __IRVING_KIRKWOOD_H__
#define __IRVING_KIRKWOOD_H__
#include "Observables.h"
#include <Aboria.h>

class PYBIND11_EXPORT VirialStress
{
    public:
        VirialStress(int _dim) : dim(_dim){};
        ~VirialStress(){};
        virtual void compute(   AboriaParticles::value_type particle_i,
                                Eigen::Ref< MatrixXd > val)
            {
                throw std::runtime_error("[ERROR] Virial Stress needs information on the j-th particle and the potential force");
            }
        virtual void compute(   AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j,
                                Eigen::Ref< MatrixXd > val)
            {
                throw std::runtime_error("[ERROR] Virial Stress needs information on the potential force");
            }
        
        virtual void compute(   AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j, 
                                const std::shared_ptr<PairPotential>& potential,
                                Eigen::Ref< MatrixXd > val)
            {
                
                //Compute pair-force!
                double forceval = potential->getPairForce();
                Eigen::Vector3d rij = potential->getRij();
                Eigen::Vector3d F =  (forceval)*(rij);
                
                //Compute virial stress
                val(0,0) += F[0]*(rij[0]);
                val(1,1) += F[1]*(rij[1]);
                val(0,1) += F[0]*(rij[1]);
                val(1,0) += F[1]*(rij[0]);
                
                if (dim > 2)
                {
                    val(2,2) += F[2]*(rij[2]);
                    val(0,2) += F[0]*(rij[2]);
                    val(2,0) += F[2]*(rij[0]);
                    val(1,2) += F[1]*(rij[2]);
                    val(2,1) += F[2]*(rij[1]);
                }
            }
        int dim;
};

class PYBIND11_EXPORT KineticStress
{
    public:
        KineticStress(int _dim) : dim(_dim){};
        ~KineticStress(){};
        
        virtual void compute(   AboriaParticles::value_type particle_i,
                        Eigen::Ref< Eigen::MatrixXd > val) 
            {
                val(0,0) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[0])*(abr::get<velocity>(particle_i)[0]);
                val(1,1) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[1])*(abr::get<velocity>(particle_i)[1]);
                val(0,1) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[0])*(abr::get<velocity>(particle_i)[1]);
                val(1,0) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[1])*(abr::get<velocity>(particle_i)[0]);
                if (dim > 2)
                {
                    val(2,2) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[2])*(abr::get<velocity>(particle_i)[2]);
                    val(0,2) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[0])*(abr::get<velocity>(particle_i)[2]);
                    val(2,0) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[2])*(abr::get<velocity>(particle_i)[0]);
                    val(1,2) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[1])*(abr::get<velocity>(particle_i)[2]);
                    val(2,1) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[2])*(abr::get<velocity>(particle_i)[1]);
                }
            }
        virtual void compute(   AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j,
                                Eigen::Ref< MatrixXd > val)
            {
                throw std::runtime_error("[ERROR] Kinetic stress only needs information on the i-th particle");
            }
        
        virtual void compute(   AboriaParticles::value_type particle_i, 
                                AboriaParticles::value_type particle_j, 
                                const std::shared_ptr<PairPotential>& potential,
                                Eigen::Ref< MatrixXd > val)
            {
                throw std::runtime_error("[ERROR] Kinetic stress only needs information on the i-th particle");
                
            }
        int dim;
};

#endif
