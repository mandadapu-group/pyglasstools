#ifndef __THERMO_PROPERTIES_H__
#define __THERMO_PROPERTIES_H__

#include <pyglasstools/MathAndTypes.h>
#include <pyglasstools/potential/PairPotential.h>
#include <Eigen/CXX11/Tensor>
namespace abr = Aboria;

class VirialStress
{
    public:
        VirialStress(){};
        ~VirialStress(){};
        
        virtual void compute(   const AboriaParticles::value_type& p_i, 
                                const AboriaParticles::value_type& p_j,
                                const Eigen::Vector3d& rij, 
                                const std::shared_ptr<PairPotential>& potential,
                                Eigen::Tensor<double, 2>& val)
        {
            
            double di = abr::get<diameter>(p_i);
            double dj = abr::get<diameter>(p_j);

            //Compute pair-force!
            double forceval = potential->getPairForce(rij, di, dj);
            Eigen::Vector3d F =  (forceval)*(rij);
            
            //Compute virial stress going through component-by-component
            for(int j = 0; j < val.dimension(1); ++j)
            {
                for (int i = 0; i < val.dimension(0); ++i)
                {
                    val(i,j) += F[i]*rij[j];
                }
            }
        }
        virtual void compute_cg(const AboriaParticles::value_type& p_i, 
                                const AboriaParticles::value_type& p_j,
                                Eigen::Vector3d rij, 
                                const std::shared_ptr<PairPotential>& potential,
                                Eigen::Tensor<double, 2>& val, double bondval)
        {
            Eigen::Tensor<double, 2> zerotensor(val);
            zerotensor.setZero();
            compute(p_i,p_j,rij,potential,zerotensor);
            val += 0.5*bondval*zerotensor;
        }
        
};

class ElasticVirialStress : public VirialStress
{
    public:
        ElasticVirialStress(){};
        ~ElasticVirialStress(){};
        
        void compute(   const AboriaParticles::value_type& p_i, 
                        const AboriaParticles::value_type& p_j,
                        const Eigen::Vector3d& rij, 
                        const std::shared_ptr<PairPotential>& potential,
                        Eigen::Tensor<double, 2>& val)
        {
            
            double di = abr::get<diameter>(p_i);
            double dj = abr::get<diameter>(p_j);

            //Compute pair-force, but now based on relative displacement!
            Eigen::Vector3d uij;
            uij  << abr::get<displacement>(p_j)[0]-abr::get<displacement>(p_i)[0],
                       abr::get<displacement>(p_j)[1]-abr::get<displacement>(p_i)[1],
                       abr::get<displacement>(p_j)[2]-abr::get<displacement>(p_i)[2];
            //py::print(newrij,abr::get<displacement>(p_j)[0],abr::get<displacement>(p_i)[0],abr::get<displacement>(p_j)[1],abr::get<displacement>(p_i)[1],
                       //abr::get<displacement>(p_j)[2],abr::get<displacement>(p_i)[2]);
            double forceval = potential->getPairForce(rij, di, dj);
            double forceval1 = potential->getPairForce(uij+rij, di, dj);
            Eigen::Vector3d F =  forceval*rij;//(urij);
            Eigen::Vector3d F1 =  forceval1*(uij+rij);//(urij);
            
            //Compute virial stress going through component-by-component
            for(int j = 0; j < val.dimension(1); ++j)
            {
                for (int i = 0; i < val.dimension(0); ++i)
                {
                    val(i,j) += F1[i]*(uij[j]+rij[j])-F[i]*rij[j];
                }
            }
        }
};

class KineticStress
{
    public:
        KineticStress(){};
        ~KineticStress(){};
        
        virtual void compute(   const AboriaParticles::value_type& p_i, 
                                Eigen::Tensor<double, 2>& val)
        {
            //Compute virial stress going through component-by-component
            for(int j = 0; j < val.dimension(1); ++j)
            {
                for (int i = 0; i < val.dimension(0); ++i)
                {
                    val(i,j) += abr::get<mass>(p_i)*(abr::get<velocity>(p_i)[i])*(abr::get<velocity>(p_i)[j]);
                }
            }
        }
        virtual void compute_cg(const AboriaParticles::value_type& p_i, 
                                Eigen::Tensor<double, 2>& val, double cgval)
        {
            compute(p_i,val);
            val = val*cgval;
        }
};

class BornStiffnessTensor
{
    public:
        BornStiffnessTensor(){};
        ~BornStiffnessTensor(){};
        
        virtual void compute(   const AboriaParticles::value_type& p_i, 
                                const AboriaParticles::value_type& p_j,
                                Eigen::Vector3d rij, 
                                const std::shared_ptr<PairPotential>& potential,
                                Eigen::Tensor<double, 4>& val)
        {
            
            //Compute pair-force!
            double di = abr::get<diameter>(p_i);
            double dj = abr::get<diameter>(p_j);
            double forceval = potential->getPairForce(rij, di, dj);
            double stiffnessval = potential->getBondStiffness(rij, di, dj);
            Eigen::Vector3d F =  (forceval)*(rij);
            Eigen::Vector3d nij = rij.normalized();

            for(int l = 0; l < val.dimension(3); ++l)
            {
                for(int k = 0; k < val.dimension(2); ++k)
                {
                    for(int j = 0; j < val.dimension(1); ++j)
                    {
                        for (int i = 0; i < val.dimension(0); ++i)
                        {
                            val(i,j,k,l) += (stiffnessval*(rij[i]) + F[i])*rij[j]*nij[k]*nij[l];
                        }
                    }
                }
            }

        }

        virtual void compute_cg(const AboriaParticles::value_type& p_i, 
                                const AboriaParticles::value_type& p_j,
                                Eigen::Vector3d rij, 
                                const std::shared_ptr<PairPotential>& potential,
                                Eigen::Tensor<double, 4>& val, double bondval)
        {
            compute(p_i,p_j,rij,potential,val);
            val = 0.5*bondval*val;
        }
};

/*
class Density
{
    public:
        Density(int _dim){};
        ~Density(){};
        virtual Eigen::MatrixXd compute(const AboriaParticles::value_type& p_i, 
                                        const AboriaParticles::value_type& p_j,
                                        Eigen::Vector3d rij, 
                                        const std::shared_ptr<PairPotential>& potential)
            {
                return Eigen::MatrixXd::Constant(1,1,1);
            }
};

*/
#endif
