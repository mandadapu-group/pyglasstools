#include "PairPotential.h"
#include "LennardJones.h" 
#include "Polydisperse12.h" 
#include "Polydisperse18.h" 
#include "PolydisperseLJ.h" 
#include "Polydisperse10.h"
#include "Polydisperse106.h"
#include "PolydisperseYukawa.h"

#include "../extern/pybind11/include/pybind11/pybind11.h"

//Typedefs several coarsegrain functions  here
typedef ShortRangePairPotential<LennardJones> PairPotentialLJ;
typedef ShortRangePairPotential<ForceShiftedLennardJones> PairPotentialForceShiftedLJ;
typedef ShortRangePairPotential<Polydisperse12> PairPotentialPoly12;
typedef ShortRangePairPotential<Polydisperse18> PairPotentialPoly18;
typedef ShortRangePairPotential<Polydisperse10> PairPotentialPoly10;
typedef ShortRangePairPotential<Polydisperse106> PairPotentialPoly106;
typedef ShortRangePairPotential<PolydisperseYukawa> PairPotentialPolyYukawa;
typedef ShortRangePairPotential<PolydisperseLJ> PairPotentialPolyLJ;


PYBIND11_MODULE(_potential, m)
{
    export_PairPotential(m);
    
    export_ShortRangePairPotential<PairPotentialLJ>(m, "PairPotentialLJ");
    export_ShortRangePairPotential<PairPotentialForceShiftedLJ>(m, "PairPotentialForceShiftedLJ");
    export_ShortRangePairPotential<PairPotentialPoly12>(m, "PairPotentialPoly12");
    export_ShortRangePairPotential<PairPotentialPoly18>(m, "PairPotentialPoly18");
    export_ShortRangePairPotential<PairPotentialPoly10>(m, "PairPotentialPoly10");
    export_ShortRangePairPotential<PairPotentialPoly106>(m, "PairPotentialPoly106");
    export_ShortRangePairPotential<PairPotentialPolyYukawa>(m, "PairPotentialPolyYukawa");
    export_ShortRangePairPotential<PairPotentialPolyLJ>(m, "PairPotentialPolyLJ");
}
