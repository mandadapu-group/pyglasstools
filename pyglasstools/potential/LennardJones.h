//(1): 1-parameter LennardJones
#ifndef __LENNARD_JONES_H__
#define __LENNARD_JONES_H__
#include <vector>

class LennardJones
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        LennardJones(double _rsq, double _rcutsq, std::vector<double> _params)
            : rsq(_rsq), rcutsq(_rcutsq), eps(_params[0])
            {
            }
        ~LennardJones(){}; 
        //! Evaluate the force and energy
        virtual double computeForce(double d_i, double d_j)
        {
            double sigma = 0.5*(d_i+d_j);
            double actualrcutsq = rcutsq*sigma*sigma;
            
            // compute the force divided by r in force_divr
            if (rsq < actualrcutsq && eps != 0)
                {
                double r2inv = 1.0/rsq;
                r2inv *= sigma*sigma;
                double r6inv = r2inv * r2inv * r2inv;
                double force_divr = 4*eps*r2inv * r6inv * (12.0*r6inv - 6.0);
                return force_divr;
                }
            else
                return 0.0;
        }
        virtual double computeRcut(double d_i, double d_j)
        {
            double sigma = 0.5*(d_i+d_j);
            return rcutsq*sigma*sigma;
        }

    protected:
        double rsq;     //!< Stored rsq from the constructor
        double rcutsq;  //!< Stored rcutsq from the constructor
        double eps;     //!< epsilon_parameter
};

//(2): 1-parameter Force-Shifted LennardJones
class ForceShiftedLennardJones : public LennardJones
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        ForceShiftedLennardJones(double _rsq, double _rcutsq, std::vector<double> _params)
            : LennardJones(_rsq, _rcutsq, _params)
            {
            }
        ~ForceShiftedLennardJones(){};        
        //! Evaluate the force and energy
        virtual double computeForce(double d_i, double d_j)
        {
            double sigma = 0.5*(d_i+d_j);
            double actualrcutsq = rcutsq*sigma*sigma;
            
            // compute the force divided by r in force_divr
            if (rsq < actualrcutsq && eps != 0)
            {
                double r2inv = 1.0/rsq;
                r2inv *= sigma*sigma;
                double r6inv = r2inv * r2inv * r2inv;
                double force_divr = 4*eps*r2inv * r6inv * (12.0*r6inv - 6.0);
                r2inv = 1.0/rcutsq;
                r6inv = r2inv * r2inv * r2inv;
                force_divr -= 4*eps*r2inv * r6inv * (12.0*r6inv - 6.0);
                return force_divr;
            }
            else
                return 0.0;
        }
};

#endif
