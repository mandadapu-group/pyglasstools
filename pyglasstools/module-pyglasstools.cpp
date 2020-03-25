#include "SimBox.h"
#include "SystemData.h"
#include "PairPotential.h"
#include "Quadrature.h"
#include "CoarseGrainFunction.h"

#include "extern/pybind11/include/pybind11/pybind11.h"

//Typedefs several pair potentials here
typedef PairPotential<LennardJones> PairPotentialLJ;
typedef PairPotential<ForceShiftedLennardJones> PairPotentialForceShiftedLJ;

//Typedefs several coarsegrain functions  here
typedef CoarseGrainFunction<Octic> OcticFunc;

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the pyglasstools python module and define the exports here.
 */

PYBIND11_MODULE(_pyglasstools, m)
{
    export_SimBox(m);
    export_SystemData(m);
    export_PairPotential<PairPotentialLJ>(m, "PairPotentialLJ");
    export_PairPotential<PairPotentialForceShiftedLJ>(m, "PairPotentialForceShiftedLJ");
    
    export_CoarseGrainFunction<OcticFunc>(m, "OcticFunc");
}
