#ifndef __POLYDISPERSE_YUKAWA_H__
#define __POLYDISPERSE_YUKAWA_H__
#include <vector>

class PolydisperseYukawa
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        PolydisperseYukawa(const double& _rsq, const double& _rcutsq, const std::vector<double>& _params)
            : rsq(_rsq), rcutsq(_rcutsq), v0(_params[0]), kappa(_params[1]), eps(_params[2]), c0(_params[3]), c1(_params[4]), c2(_params[5])
            {                
            }
        ~PolydisperseYukawa(){}; 
        //! Evaluate energy
        virtual double computeEnergy(const double& d_i, const double& d_j) const
        {
            double sigma = 0.5*(d_i+d_j)*(1-eps*fabs(d_i-d_j));
            double actualrcutsq = rcutsq*sigma*sigma;
            
            // compute the force divided by r in force_divr
            if (rsq < actualrcutsq)
                {
                    double r2inv = 1.0/rsq;
                    r2inv *= sigma*sigma;
                    double _rsq = rsq/(sigma*sigma);
                    double rinv = sqrt(r2inv);
                    double rroot = sqrt(_rsq);
                    return v0*exp(-rroot*kappa)*rinv+c0+c1*_rsq + c2*_rsq*_rsq;
                }
            else
                return 0.0;
        }
        
        //! Evaluate the force and energy
        virtual double computeForce(const double& d_i, const double& d_j) const
        {
            double sigma = 0.5*(d_i+d_j)*(1-eps*fabs(d_i-d_j));
            double actualrcutsq = rcutsq*sigma*sigma;
            
            // compute the force divided by r in force_divr
            if (rsq < actualrcutsq)
                {
                    double r2inv = 1.0/rsq;
                    r2inv *= sigma*sigma;
                    double _rsq = rsq/(sigma*sigma);
                    double rinv = sqrt(r2inv);
                    double rroot = sqrt(_rsq);
                    return ((exp(-kappa*rroot)*(rinv + kappa)*v0)*r2inv - 2*c1 - 4*c2*_rsq)/(sigma*sigma);
                }
            else
                return 0.0;
        }
        
        //! Evaluate the force and energy
        virtual double computeStiffness(const double& d_i, const double& d_j) const
        {
            double sigma = 0.5*(d_i+d_j)*(1-eps*fabs(d_i-d_j));
            double actualrcutsq = rcutsq*sigma*sigma;
            
            // compute the secon derivative of pair potential
            if (rsq < actualrcutsq)
                {
                    double r2inv = 1.0/rsq;
                    r2inv *= sigma*sigma;
                    double _rsq = rsq/(sigma*sigma);
                    double rinv = sqrt(r2inv);
                    double rroot = sqrt(_rsq);
                    return (v0*exp(-kappa*rroot)*(2*r2inv + 2*kappa*rinv + kappa*kappa)*rinv + 2*c1 + 12*c2*_rsq)/(sigma*sigma);
                }
            else
                return 0.0;
        }
        
        virtual double computeRcut(const double& d_i, const double& d_j) const
        {
            double sigma = 0.5*(d_i+d_j)*(1-eps*fabs(d_i-d_j));
            return rcutsq*sigma*sigma;
        }

    protected:
        double rsq;     //!< Stored rsq from the constructor
        double rcutsq;  //!< Stored rcutsq from the constructor
        double v0;
        double kappa;
        double eps;     //!< epsilon_parameter
        double c0;
        double c1;
        double c2;
};

#endif
