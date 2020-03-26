#include "SimBox.h"
#include "ParticleSystem.h"
#include "PairPotential.h"
#include "Quadrature.h"
#include "CoarseGrainFunction.h"

#include "extern/pybind11/include/pybind11/pybind11.h"

//Typedefs several pair potentials here
typedef ShortRangePairPotential<LennardJones> PairPotentialLJ;
//typedef ParticleSystem<PairPotentialLJ> LennardJonesSystem;
typedef ShortRangePairPotential<ForceShiftedLennardJones> PairPotentialForceShiftedLJ;

//Typedefs several coarsegrain functions  here
typedef CoarseGrainFunction<Octic> OcticFunc;

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the pyglasstools python module and define the exports here.
 */

PYBIND11_MODULE(_pyglasstools, m)
{
    export_SimBox(m);
    export_PairPotential(m);
    export_ParticleSystem<ParticleSystem>(m, "ParticleSystem");
    
    export_ShortRangePairPotential<PairPotentialLJ>(m, "PairPotentialLJ");
    export_ShortRangePairPotential<PairPotentialForceShiftedLJ>(m, "PairPotentialForceShiftedLJ");
    
    
    export_CoarseGrainFunction<OcticFunc>(m, "OcticFunc");
}
