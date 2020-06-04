//#include "SpectraHessian.h"
#include "SLEPcHessian.h"
#include <pybind11/pybind11.h>

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the observables python module and define the exports here.
 */

typedef PETScVectorField<2> PETScVectorField2D;
typedef PETScVectorField<3> PETScVectorField3D;
typedef PETScGlobalProperty<1,2> GlobalVector2D;   
typedef PETScGlobalProperty<1,3> GlobalVector3D;   

PYBIND11_MODULE(_nonaffine, m)
{
    int external_init = slepc::initialize();

    if (!external_init)
    {
        Py_AtExit(slepc::finalize);
    }
    export_PETScVectorFieldBase(m);
    export_PETScGlobalPropertyBase(m);
    export_PETScGlobalProperty< GlobalVector2D >(m,"GlobalVector2D");    
    export_PETScGlobalProperty< GlobalVector3D >(m,"GlobalVector3D");    
    export_PETScVectorField< PETScVectorField2D >(m,"PETScVectorField2D");
    export_PETScVectorField< PETScVectorField3D >(m,"PETScVectorField3D");
    export_PETScManager(m);
    export_HessianManager(m);
    export_SLEPcHessian(m);
    //export_SpectraHessian(m);
}
