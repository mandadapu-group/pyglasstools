#include "PETScManager.h"
#include "PETScHessian.h"
#include "PETScVectorField.h"
#include "PETScCalculatorBase.h"
#include "PETScLinearResponse.h"
#include "SLEPcNMA.h"

#include <slepceps.h>
#include <pybind11/pybind11.h>

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the observables python module and define the exports here.
 */

namespace slepc
{
    PetscErrorCode initialize()
    {
        PetscErrorCode ierr = SlepcInitialize(NULL,NULL,(char*)0,NULL);
        return ierr;
    }

    void finalize()
    {
        SlepcFinalize();
    }
};

typedef PETScVectorField<2> PETScVectorField2D;
typedef PETScVectorField<3> PETScVectorField3D;

PYBIND11_MODULE(_elasticitypetsc, m)
{
    int external_init = slepc::initialize();

    if (!external_init)
    {
        Py_AtExit(slepc::finalize);
    }
    //Export PETSc Managers
    export_PETScManager(m);
    
    //Export PETScObservables
    export_PETScVectorFieldBase(m);
    export_PETScVectorField< PETScVectorField2D >(m,"PETScVectorField2D");
    export_PETScVectorField< PETScVectorField3D >(m,"PETScVectorField3D");
    
    //Export PETSc Hessian Objects
    export_PETScHessianBase(m);
    export_PETScHessian(m);
    
    //Export PETSc Calculators
    export_PETScCalculatorBase(m);
    export_PETScLinearResponse(m);
    export_SLEPcNMA(m);
}
