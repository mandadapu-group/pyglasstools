#ifndef __MODEL_LJ_H__
#define __MODEL_LJ_H__

class ModelLJ
    {
    public:
        //! Define the parameter type used by this pair potential evaluator
        typedef double param_type[2];

        //! Constructs the pair potential evaluator
        /*! \param _rsq Squared distance between the particles
            \param _rcutsq Squared distance at which the potential goes to 0
            \param _params Per type pair parameters of this potential
        */
        DEVICE EvaluatorPairLJ(double _rsq, double _rcutsq, const param_type& _params)
            : rsq(_rsq), rcutsq(_rcutsq), lj1(_params[0]), lj2(_params[1])
            {
            }

        //! LJ doesn't use diameter
        DEVICE static bool needsDiameter() { return false; }
        //! Accept the optional diameter values
        /*! \param di Diameter of particle i
            \param dj Diameter of particle j
        */
        DEVICE void setDiameter(double di, double dj) { }

        //! LJ doesn't use charge
        DEVICE static bool needsCharge() { return false; }
        //! Accept the optional diameter values
        /*! \param qi Charge of particle i
            \param qj Charge of particle j
        */
        DEVICE void setCharge(double qi, double qj) { }

        //! Evaluate the force and energy
        /*! \param force_divr Output parameter to write the computed force divided by r.
            \param pair_eng Output parameter to write the computed pair energy
            \param energy_shift If true, the potential must be shifted so that V(r) is continuous at the cutoff
            \note There is no need to check if rsq < rcutsq in this method. Cutoff tests are performed
                  in PotentialPair.

            \return True if they are evaluated or false if they are not because we are beyond the cutoff
        */
        double computeforce()
            {
            // compute the force divided by r in force_divr
            if (rsq < rcutsq && lj1 != 0)
                {
                double r2inv = 1.0/rsq;
                double r6inv = r2inv * r2inv * r2inv;
                force_divr= r2inv * r6inv * (12.0*lj1*r6inv - 6.0*lj2);

                return true;
                }
            else
                return 0.0;
            }

    protected:
        double rsq;     //!< Stored rsq from the constructor
        double rcutsq;  //!< Stored rcutsq from the constructor
        double lj1;     //!< lj1 parameter extracted from the params passed to the constructor
        double lj2;     //!< lj2 parameter extracted from the params passed to the constructor
    };


#endif // __PAIR_EVALUATOR_LJ_H__
