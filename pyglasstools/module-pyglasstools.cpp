#include "SimBox.h"
#include "AllPluginPairPotentials.h"
#include "PairPotential.h"

#include "extern/pybind11/include/pybind11/pybind11.h"
//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the hoomd python module and define the exports here.
 */

PYBIND11_MODULE(_pyglasstools, m)
    {
    export_SimBox(m);
    export_PotentialPair<PotentialPairLJPlugin>(m, "PotentialPairLJPlugin");
    }
