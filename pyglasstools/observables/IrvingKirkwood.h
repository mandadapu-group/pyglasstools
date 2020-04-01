#ifndef __IRVING_KIRKWOOD_H__
#define __IRVING_KIRKWOOD_H__

#include <pyglasstools/MathAndTypes.h>
namespace abr = Aboria;

class VirialStress
{
    public:
        VirialStress(int _dim) : dim(_dim){};
        ~VirialStress(){};
        virtual Eigen::MatrixXd compute(   const AboriaParticles::value_type& particle_i)
            {
                throw std::runtime_error("[ERROR] Virial Stress needs information on the j-th particle and the potential force");
                return Eigen::MatrixXd::Zero(dim,dim);
            }
        virtual Eigen::MatrixXd compute(   const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j)
            {
                throw std::runtime_error("[ERROR] Virial Stress needs information on the potential force");
                return Eigen::MatrixXd::Zero(dim,dim);
            }
        
        virtual Eigen::MatrixXd compute(   const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j, 
                                const std::shared_ptr<PairPotential>& potential)
            {
                
                //Compute pair-force!
                double forceval = potential->getPairForce();
                Eigen::Vector3d rij = potential->getRij();
                Eigen::Vector3d F =  (forceval)*(rij);
                Eigen::MatrixXd val =  Eigen::MatrixXd::Zero(dim,dim);
                
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
                return val;
            }
        int dim;
};

class KineticStress
{
    public:
        KineticStress(int _dim) : dim(_dim){};
        ~KineticStress(){};
        
        virtual Eigen::MatrixXd compute(   const AboriaParticles::value_type& particle_i)
            {
                Eigen::MatrixXd val =  Eigen::MatrixXd::Zero(dim,dim);
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
                return val;
            }
        virtual Eigen::MatrixXd compute(   const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j)
            {
                throw std::runtime_error("[ERROR] Kinetic stress only needs information on the i-th particle");
                return Eigen::MatrixXd::Zero(dim,dim);
            }
        
        virtual Eigen::MatrixXd compute(   const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j, 
                                const std::shared_ptr<PairPotential>& potential)
            {
                throw std::runtime_error("[ERROR] Kinetic stress only needs information on the i-th particle");
                return Eigen::MatrixXd::Zero(dim,dim);
                
            }
        int dim;
};

class Density
{
    public:
        Density(int _dim){};
        ~Density(){};
        
        virtual Eigen::MatrixXd compute(   const AboriaParticles::value_type& particle_i)
            {
                return Eigen::MatrixXd::Constant(1,1,1);
            }
        virtual Eigen::MatrixXd compute(   const AboriaParticles::value_type& particle_i, 
                                        const AboriaParticles::value_type& particle_j)
            {
                throw std::runtime_error("[ERROR] Density only needs information on the i-th particle");
                return Eigen::MatrixXd::Zero(1,1);
            }
        
        virtual Eigen::MatrixXd compute(   const AboriaParticles::value_type& particle_i, 
                                            const AboriaParticles::value_type& particle_j, 
                                            const std::shared_ptr<PairPotential>& potential)
            {
                throw std::runtime_error("[ERROR] Density only needs information on the i-th particle");
                return Eigen::MatrixXd::Zero(1,1);
            }
};

#endif
