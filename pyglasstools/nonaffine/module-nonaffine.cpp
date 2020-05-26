#include "SpectraHessian.h"
#include "SLEPcHessian.h"
#include <pybind11/pybind11.h>

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the observables python module and define the exports here.
 */

PYBIND11_MODULE(_nonaffine, m)
{
    int external_init = slepc::initialize();

    if (!external_init)
    {
        Py_AtExit(slepc::finalize);
    }
    
    export_PETScManager(m);
    export_HessianManager(m);
    export_SLEPcHessian(m);
    export_SpectraHessian(m);
}
