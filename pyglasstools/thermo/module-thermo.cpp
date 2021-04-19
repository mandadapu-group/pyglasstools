#include "ThermoProperty.h"
#include "ThermoCalculator.h"
#include <pyglasstools/ListofCommonObservables.h>
#include <pybind11/pybind11.h>

//Typedefs several observables here
typedef ForceScalarProperty< 2, PotentialEnergy > GlobalPotentialEnergy2D;
typedef ForceScalarProperty< 3, PotentialEnergy > GlobalPotentialEnergy3D;

typedef ForceProperty<2, 2, VirialStress > GlobalVirialStress2D;
typedef ForceProperty<2, 3, VirialStress > GlobalVirialStress3D;

typedef ForceProperty<2, 2, ElasticVirialStress > GlobalElasticVirialStress2D;
typedef ForceProperty<2, 3, ElasticVirialStress > GlobalElasticVirialStress3D;

typedef LocalProperty<2, 2, KineticStress > GlobalKineticStress2D;
typedef LocalProperty<2, 3, KineticStress > GlobalKineticStress3D;

typedef ForceProperty<4, 2, BornStiffnessTensor > GlobalBornTensor2D;
typedef ForceProperty<4, 3, BornStiffnessTensor > GlobalBornTensor3D;

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the observables python module and define the exports here.
 */

PYBIND11_MODULE(_thermo, m)
{
    export_ThermoCalculator(m);
    export_ThermoPropertyBase(m);
    
    export_ThermoProperty< GlobalPotentialEnergy2D >(m, "GlobalPotentialEnergy2D");
    export_ThermoProperty< GlobalPotentialEnergy3D >(m, "GlobalPotentialEnergy3D");
    
    export_ThermoProperty< GlobalVirialStress2D >(m, "GlobalVirialStress2D");
    export_ThermoProperty< GlobalVirialStress3D >(m, "GlobalVirialStress3D");
    
    export_ThermoProperty< GlobalElasticVirialStress2D >(m, "GlobalElasticVirialStress2D");
    export_ThermoProperty< GlobalElasticVirialStress3D >(m, "GlobalElasticVirialStress3D");
    
    export_ThermoProperty< GlobalBornTensor2D >(m, "GlobalBornTensor2D");
    export_ThermoProperty< GlobalBornTensor3D >(m, "GlobalBornTensor3D");
    
    export_ThermoProperty< GlobalKineticStress2D >(m, "GlobalKineticStress2D");
    export_ThermoProperty< GlobalKineticStress3D >(m, "GlobalKineticStress3D");
}
