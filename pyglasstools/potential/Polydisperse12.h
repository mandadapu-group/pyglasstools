#ifndef __POLYDISPERSE_12_H__
#define __POLYDISPERSE_12_H__
#include <vector>

class Polydisperse12
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        Polydisperse12(double _rsq, double _rcutsq, std::vector<double> _params)
            : rsq(_rsq), rcutsq(_rcutsq), v0(_params[0]), eps(_params[1])
            {
                c0 =  -28.0*v0/pow(_rcutsq,6);
                c1 =  48.0*v0/pow(_rcutsq,7);
                c2 =  -21.0*v0/pow(_rcutsq,8);
            }
        ~Polydisperse12(){}; 
        //! Evaluate the force and energy
        virtual double computeForce(double d_i, double d_j)
        {
            double sigma = 0.5*(d_i+d_j)*(1-eps*fabs(d_i-d_j));
            double actualrcutsq = rcutsq*sigma*sigma;
            
            // compute the force divided by r in force_divr
            if (rsq < actualrcutsq)
                {
                    double r2inv = 1.0/rsq;
                    r2inv *= sigma*sigma;
                    double _rsq = rsq/(sigma*sigma);
                    double r6inv = r2inv * r2inv * r2inv;
                    double force_divr = 12.0*v0*r2inv*r6inv*r6inv-2.0*c1 -4.0*c2*_rsq;
                    return force_divr;
                }
            else
                return 0.0;
        }
        virtual double computeRcut(double d_i, double d_j)
        {
            double sigma = 0.5*(d_i+d_j)*(1-eps*fabs(d_i-d_j));
            return rcutsq*sigma*sigma;
        }

    protected:
        double rsq;     //!< Stored rsq from the constructor
        double rcutsq;  //!< Stored rcutsq from the constructor
        double v0;
        double eps;     //!< epsilon_parameter
        double c0;
        double c1;
        double c2;
};

#endif
