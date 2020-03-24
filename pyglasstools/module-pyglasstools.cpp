#include "SimBox.h"
#include "PairPotential.h"

#include "extern/pybind11/include/pybind11/pybind11.h"

typedef PairPotential<LennardJones> PairPotentialLJ;
typedef PairPotential<ForceShiftedLennardJones> PairPotentialForceShiftedLJ;
//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the pyglasstools python module and define the exports here.
 */

PYBIND11_MODULE(_pyglasstools, m)
{
    export_SimBox(m);
    export_PairPotential<PairPotentialLJ>(m, "PairPotentialLJ");
    export_PairPotential<PairPotentialForceShiftedLJ>(m, "PairPotentialForceShiftedLJ");
}
