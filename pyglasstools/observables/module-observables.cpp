#include "Observables.h"
#include "IrvingKirkwood.h"

#include "../extern/pybind11/include/pybind11/pybind11.h"

//Typedefs several observables here
typedef GlobalObservable<VirialStress> GlobalVirialStress;
typedef GlobalObservable<KineticStress> GlobalKineticStress;
typedef LocalObservable<VirialStress> VirialStressField;
typedef LocalObservable<KineticStress> KineticStressField;
typedef LocalObservable<Density> DensityField;
//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the observables python module and define the exports here.
 */

PYBIND11_MODULE(_observables, m)
{
    export_Observable(m);
    export_GlobalObservable<GlobalVirialStress>(m, "GlobalVirialStress");
    export_GlobalObservable<GlobalKineticStress>(m, "GlobalKineticStress");
    
    export_LocalObservable<VirialStressField>(m, "VirialStressField");
    export_LocalObservable<KineticStressField>(m, "KineticStressField");
    export_LocalObservable<DensityField>(m, "DensityField");
}
