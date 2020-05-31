#include <cmath>
class Rect
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        Rect(double _rcut)
            : rcut(_rcut)
            {}
        ~Rect(){}; 
        virtual void setRcut(double _rcut){rcut=_rcut;}; 
        
        //! Evaluate the force and energy
        virtual double compute(double dr_sq)
        {
            double rcutsq = rcut*rcut;
            
            // compute the force divided by r in force_divr
            if (dr_sq < rcutsq)
            {
                return 1/(M_PI*rcutsq);
            }
            else
                return 0.0;
        }

    protected:
        double rcut;  //!< Stored rcutsq from the constructor
};
