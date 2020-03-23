#include <vector>
#include <array>
#include <math>
#include "PairPotential.h"

class CoarseGrainFunction
{
    public:
        double compute(std:array<double,3> x, std::array<double,3> rprime)
        dr = PBC(x-rprime,boxsize)#/rcut
        absr = np.sqrt(np.sum(dr**2))/rcut#-rcut
        if absr < 1:
            return 15/(8*np.pi*rcut**2)*(1-2*absr**4+absr**8)
        else:
            return 0
    private:
        double boxsize;
        double cg_rcut;
        
    //@nb.njit(nogil=True,fastmath=True,cache=True)
    def Delta(x,rprime,rcut,boxsize):
}

class BondFunction
{
    @nb.njit(nogil=True,fastmath=True,cache=True)
    def bond(s,x,ri,dr,rcut,boxsize):
        return Delta(x,ri+s*dr,rcut,boxsize)

    weightquad = np.array([0.568888889,0.47862867,0.47862867,0.236926885,0.236926885])
    xquad = np.array([0,-0.53846931,0.53846931,-0.906179846, 0.906179846])     
    @nb.njit(nogil=True,fastmath=True,cache=True)
    def integrate_bond(x,ri,dr,rcut,boxsize):
        y = (xquad+1)/2.0 
        val = 0
        for i in range(5):
            val += 0.5*(weightquad[i]*bond(y[i],x,ri,dr,rcut,boxsize))
        return val
}

class IKField
{
    public:
        CoarseGrainFields();
        ~CoarseGrainFields(){};
        void computeCGfield();
        double applyPBC(std::array<double, 3> dr);
    private:
        //All Possible Coarse-Grained Fields
        std::vector<double> density;
        std::vector< std::array<double, 6> > virial_stress;
        std::vector< std::array<double, 6> > kinetic_stress;
        
        //Atomic properties, fed in to the system
        std::vector<double> diameter;
        std::vector< std::array<double, 3> > m_gridposition;
        std::vector< std::array<double, 3> > m_atomposition;
        std::vector< std::array<double, 3> > m_atomvelocity;
        std::vector< std::vector<unsigned int> > m_coarsegrain_pair;  
        std::vector< std::vector<unsigned int> > m_force_pair;  

        std::shared_ptr< PairPotential > m_potential;     //!< Patchy Interaction
        std::shared_ptr< BondFunction > m_bondfunc; 
        std::shared_ptr< CoarseGrainFunction > m_cgfunc; 
        
        double m_boxsize;
        double m_force_rcut;
        double m_cg_rcut;
        int m_dim;
        bool m_hasvelocities;
}

double applyPBC(std::array<double, 3> dr)
{
    for i in range(dx.shape[0]):
    {
        if (i > m_dim)
        {
            continue;
        }
        else
        {
            if (dr[i] > m_boxsize*0.5)
                dr[i] -= m_boxsize;
            else if (dx[i] <= -m_boxsize*0.5):
                dr[i] += m_boxsize;
            else
                continue;
        }
    }
    return dr
}

void CoarseGrainFields::computeCGField()
{
    for(unsigned int n = 0, n < m_atomposition.size(), n++)
    {
        std::vector<int> pair_ni(m_coarsegrain_pair[n].begin(),m_coarsegrain_pair[n].end());
        
        for(unsigned int i = 0, i < pair_ni.size(), i++)
        {
            unsigned int id_i = pair_ni[i];
            rho[n] += m_cgfunc->computeCGFunc(  m_gridposition[n],
                                                m_atomposition[i],
                                                m_cg_rcut,
                                                m_boxsize);
            if (m_hasvelocities)
            {
            }
            std::vector<int> pair_ij(m_force_pair[i].begin(),m_force_pair[i].end());
            
            for(unsigned int j = 0, j < pair_ij.size(), i++)
            {
                unsigned int id_j = pair_ij[j];
                if (id_i != id_j)
                {
                    std::array<double, 3> dr;
                    dr = applyPBC(m_atomposition[id_j]-m_atomposition[id_i],m_boxsize); 
                    if (m_pairpotential->checkRange(dr) )
                    {
                        std::array<double, 3> force_ij;
                        m_pairpotential->setDiameters(diameters[id_i],diameters[id_j])
                        force_ij = m_pairpotential->computeforce(); 
                        
                        if (m_bondfunc->checkRange(dr))
                        {
                            double bond_val = m_bondfunc->compute(m_gridposition[n],
                                                                  m_atomposition[i],
                                                                  dr);

                            virial_stress[n][0] += 0.5*F[0]*dr[0]*bond_val;
                            virial_stress[n][1] += 0.5*F[1]*dr[1]*bond_val;
                            virial_stress[n][2] += 0.5*F[0]*dr[1]*bond_val;
                            if (m_dim > 2)
                            {
                                virial_stress[n][3] += 0.5*F[2]*dr[2]*bond_val;
                                virial_stress[n][4] += 0.5*F[0]*dr[2]*bond_val;
                                virial_stress[n][5] += 0.5*F[1]*dr[2]*bond_val;
                            } //end of if m_dim > 2
                        } //end of m_bondfunc->checkRange(dr)
                    } //end of m_pairpotential->checkRange(dr) 
                } //end of if id_i != id_j
            } //end of for loop j
        } //end of for loop i
    } //end of for loop n
}
