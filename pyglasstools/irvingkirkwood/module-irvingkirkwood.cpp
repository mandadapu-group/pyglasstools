#include "CoarseGrainedField.h"
#include "IrvingKirkwood.h"
#include <pyglasstools/ListofCommonObservables.h>

#include <pybind11/pybind11.h>

//Typedefs several observables here
typedef ForceField< 2, 2, VirialStress> VirialStressField2D;
typedef ForceField< 2, 3, VirialStress> VirialStressField3D;
typedef LocalField< 2, 2, KineticStress> KineticStressField2D;
typedef LocalField< 2, 3, KineticStress> KineticStressField3D;
//typedef CoarseGrainedField<Density> DensityField;

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the observables python module and define the exports here.
 */

PYBIND11_MODULE(_irvingkirkwood, m)
{
    export_IrvingKirkwood(m);
    export_CoarseGrainedFieldBase(m);
    export_CoarseGrainedField< VirialStressField2D >(m, "VirialStressField2D");
    export_CoarseGrainedField< KineticStressField2D >(m, "KineticStressField2D");
    export_CoarseGrainedField< VirialStressField3D >(m, "VirialStressField3D");
    export_CoarseGrainedField< KineticStressField3D >(m, "KineticStressField3D");
    //export_CoarseGrainedField<DensityField>(m, "DensityField");
}
