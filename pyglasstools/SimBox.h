//Look at SimBox.h
#ifndef __SIMBOX_H__
#define __SIMBOX_H__

#include <array>

class SimBox
{
    public:
        //A dummy constructor, just to have something setup when virtually no argument is given
        SimBox()
        {
            m_L.fill(0);
            m_Lmin.fill(0);
            m_Lmax.fill(0);
            m_Linv.fill(0);
            m_periodic.fill(true);
        };

        //! Constructs a box from -Len/2 to Len/2
        /*! \param Len Length of one side of the box
            \post Box ranges from \c -Len/2 to \c +Len/2 in all 3 dimensions
            \post periodic = (1,1,1)
        */
        SimBox(double Len)
        {
            m_L.fill(Len);
            setL(m_L);
            m_periodic.fill(true); 
        };

        //! Constructs a box from -Len_x/2 to Len_x/2 for each dimension
        /*! \param Len_x Length of the x dimension of the box
            \param Len_y Length of the x dimension of the box
            \param Len_z Length of the x dimension of the box
            \post periodic = (1,1,1)
        */
        SimBox(double Len_x, double Len_y, double Len_z)
        {
            std::array<double, 3> inputL = {Len_x, Len_y, Len_z};
            setL(inputL);
            m_periodic.fill(1);
        };
        
        ~SimBox(){};

        //! Get the periodic flags
        /*! \return Periodic flags
        */
        std::array<bool,3> getPeriodic() const
        {
            return m_periodic;
        };

        //! Set the periodic flags
        /*! \param periodic Flags to set
            \post Period flags are set to \a periodic
            \note It is invalid to set 1 for a periodic dimension where lo != -hi. This error is not checked for.
        */
        void setPeriodic(std::array<bool, 3> periodic)
        {
            m_periodic = periodic;
        };


        //! Update the box length
        /*! \param L new box length in each direction
        */
        void setL(const std::array<double, 3>& L)
        {
            m_Lmax = L/2.0;//Scalar(2.0);
            m_Lmin = -m_Lmax;
            m_Linv = 1.0/L;
            m_L = L;
        }
        //! Get the length of the box in each direction
        /*! \returns The length of the box in each direction (hi - lo)
        */
        std::array<double,3> getL() const
        {
            return m_L;
        }

        //! Apply periodic boundary conditions to a vector
        void applyPBC(std::array<double,3>& v)
        {
            std::array<double, 3> L = getL();
        
            //A bunch of if else statements here
            //z-direction
            if (m_periodic[2])
            {
                if (v[2] >= m_Lmax[2])
                    v[2] -= L[2];
                else if (v[2] < m_Lmin[2])
                    v[2] += L[2];
            }
            //y-direction
            if (m_periodic[1])
            {
                if (v[1] >= m_Lmax[1])
                    v[1] -= L[1];
                else if (v[1] < m_Lmin[1])
                    v[1] += L[1];
            }
            //x-direction
            if (m_periodic[0])
            {
                if (v[0] >= m_Lmax[0])
                    v[0] -= L[0];
                else if (v[0] < m_Lmin[0])
                    v[0] += L[0];
            }
        }

        //! Get the volume of the box
        /*! \returns the volume
         *  \param twod If true, return the area instead of the volume
         */
        double getVolume(bool twod=false) const
        {
            if (twod)
                return m_L[0]*m_L[1];
            else
                return m_L[0]*m_L[1]*m_L[2];
        }

    private:
        std::array<double, 3> m_Lmin;       //!< L precomputed (used to avoid subtractions in boundary conditions)
        std::array<double, 3> m_Lmax;       //!< minimum value of L, per coordinate precomputed
        std::array<double, 3> m_L;       //!< L precomputed (used to avoid subtractions in boundary conditions)
        std::array<double, 3> m_Linv;       //!< 1.L precomputed (used to avoid subtractions in boundary conditions)
        std::array<bool, 3> m_periodic; //!< 0/1 in each direction to tell if the box is periodic in that direction
};

//an export function here
#endif // __SIMBOX_H__
