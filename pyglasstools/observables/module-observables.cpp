#include "Observables.h"
#include "IrvingKirkwood.h"

#include "../extern/pybind11/include/pybind11/pybind11.h"

//Typedefs several observables here
typedef GlobalObservable<VirialStress> GlobVirialStress;
typedef GlobalObservable<KineticStress> GlobKineticStress;
//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the observables python module and define the exports here.
 */

PYBIND11_MODULE(_observables, m)
{
    export_Observable(m);
    export_GlobalObservable<GlobVirialStress>(m, "GlobVirialStress");
    export_GlobalObservable<GlobKineticStress>(m, "GlobKineticStress");
}
