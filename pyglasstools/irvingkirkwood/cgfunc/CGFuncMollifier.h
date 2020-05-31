#include <cmath>

class Mollifier
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        Mollifier(double _rcut)
            : NORM_CONST(13.46842098743079272180855772032086697297031125109170011102/(2.0*M_PI*_rcut*_rcut)),
              rcut(_rcut)
            {
            }
        ~Mollifier(){}; 
        virtual void setRcut(double _rcut){rcut=_rcut;}; 

        //! Evaluate the force and energy
        virtual double compute(double dr_sq)
        {
            double rcutsq = rcut*rcut;
            
            // compute the force divided by r in force_divr
            if (dr_sq < rcutsq)
            {
                double r2 = dr_sq/rcutsq;
                return NORM_CONST*exp(-1.0/(1-r2));
            }
            else
                return 0.0;
        }

    protected:
        double NORM_CONST;
        double rcut;  //!< Stored rcutsq from the constructor
};
