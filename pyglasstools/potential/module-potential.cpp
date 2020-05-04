#include "PairPotential.h"
#include "LennardJones.h" 
#include "Polydisperse12.h" 
#include "Polydisperse18.h" 

#include "../extern/pybind11/include/pybind11/pybind11.h"

//Typedefs several coarsegrain functions  here
typedef ShortRangePairPotential<LennardJones> PairPotentialLJ;
typedef ShortRangePairPotential<ForceShiftedLennardJones> PairPotentialForceShiftedLJ;
typedef ShortRangePairPotential<Polydisperse12> PairPotentialPoly12;
typedef ShortRangePairPotential<Polydisperse18> PairPotentialPoly18;


PYBIND11_MODULE(_potential, m)
{
    export_PairPotential(m);
    
    export_ShortRangePairPotential<PairPotentialLJ>(m, "PairPotentialLJ");
    export_ShortRangePairPotential<PairPotentialForceShiftedLJ>(m, "PairPotentialForceShiftedLJ");
    export_ShortRangePairPotential<PairPotentialPoly12>(m, "PairPotentialPoly12");
    export_ShortRangePairPotential<PairPotentialPoly18>(m, "PairPotentialPoly18");
}
