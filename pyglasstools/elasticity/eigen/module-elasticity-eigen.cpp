#include "EigenManager.h"
#include "EigenHessian.h"
#include "EigenVectorField.h"
#include "EigenLinearResponse.h"

//#include "PETScForceDipole.h"
//#include "SpectraNMA.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the observables python module and define the exports here.
 */
typedef EigenVectorField<2> EigenVectorField2D;
typedef EigenVectorField<3> EigenVectorField3D;

PYBIND11_MODULE(_elasticityeigen, m)
{
    //Export PETSc Managers
    export_EigenManager(m);
        
    export_EigenVectorFieldBase(m);
    export_EigenVectorField< EigenVectorField2D >(m,"EigenVectorField2D");
    export_EigenVectorField< EigenVectorField3D >(m,"EigenVectorField3D");
    
    export_EigenHessianBase(m);
    export_EigenHessian< EigenHessian<2> >(m, "EigenHessian2D");
    export_EigenHessian< EigenHessian<3> >(m, "EigenHessian3D");
    
    //Export PETSc Calculators
    export_EigenLinearResponse(m);
}
