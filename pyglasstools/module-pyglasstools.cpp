#include "SimBox.h"
#include "ParticleSystem.h"
#include "GlobalCalculator.h"

#include "extern/pybind11/include/pybind11/pybind11.h"

//Typedefs several pair potentials here

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the pyglasstools python module and define the exports here.
 */

PYBIND11_MODULE(_pyglasstools, m)
{
    export_SimBox(m);
    
    export_ParticleSystem(m);
    
    export_GlobalCalculator(m);    
    
    export_LocalCalculator(m);    
}
