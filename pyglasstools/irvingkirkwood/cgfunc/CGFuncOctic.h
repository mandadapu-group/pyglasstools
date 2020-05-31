class Octic
{
    public:
        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        Octic(double _rcut)
            : rcut(_rcut)
            {
            }
        ~Octic(){};
        virtual void setRcut(double _rcut){rcut=_rcut;}; 
        //! Evaluate the force and energy
        virtual double compute(double dr_sq)
        {
            double rcutsq = rcut*rcut;
            
            // compute the force divided by r in force_divr
            if (dr_sq < rcutsq)
            {
                double r4 = dr_sq*dr_sq/(rcutsq*rcutsq);
                double r8 = r4*r4;
                return 15.0/(8*M_PI*rcutsq)*(1-2*r4+r8);
            }
            else
                return 0.0;
        }

    protected:
        double rcut;  //!< Stored rcutsq from the constructor
};
