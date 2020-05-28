#include "ThermoProperty.h"
#include "ThermoCalculator.h"
#include <pyglasstools/ListofCommonObservables.h>
#include <pybind11/pybind11.h>

//Typedefs several observables here
typedef ThermoProperty<VirialStress> GlobalVirialStress;
typedef ThermoProperty<KineticStress> GlobalKineticStress;
typedef ThermoProperty<Density> GlobalDensity;
typedef ThermoProperty<BornStiffnessTensor> GlobalBornTensor;

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the observables python module and define the exports here.
 */

PYBIND11_MODULE(_thermo, m)
{
    //export_Observable(m);
    export_ThermoCalculator(m);
    export_ThermoProperty<GlobalVirialStress>(m, "GlobalVirialStress");
    export_ThermoProperty<GlobalKineticStress>(m, "GlobalKineticStress");
    export_ThermoProperty<GlobalBornTensor>(m, "GlobalBornTensor");
    export_ThermoProperty<GlobalDensity>(m, "GlobalDensity");
    
}
