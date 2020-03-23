//a template class for constructing lots of pair potentials
#ifndef __PAIR_POTENTIAL_H__
#define __PAIR_POTENTIAL_H__

#include <vector>
#include <array>
#include <pyglasstools/extern/pybind/include/pybind11/pybind11.h>

// A Class that can store multiples
template < class Model >
class PYBIND11_EXPORT PairPotential
{
    public:
        typedef typename model::param_type param_type;
        
        PairPotential();
        virtual ~PairPotential(){};
        
        virtual double getPairForce();
        {
            double rsq_ij;
            rsq_ij = r_ij[0]*r_ij[0]+r_ij[1]*r_ij[1]+r_ij[2]*r_ij[2];
            double rsq_cut = rcut*rcut;
            
            Model model(rsq, rcutsq, d_i, d_j);
            return model.computeforce()
        }
        virtual void setDiameters(double _d_i, double _d_j);
        {
            d_i = _d_i;
            d_j = _d_j;
        }

    private:
        std::array<double, 3> r_ij;
        double d_i;
        double d_j; 
};

PairPotential::PairPotential
