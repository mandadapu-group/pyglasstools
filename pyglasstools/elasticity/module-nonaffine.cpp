#include "NonAffineManager.h"

//#include "PETScHessian.h"
#include "SLEPcHessian.h"
//#include "SpectraHessian.h"

#include "GlobalProperty.h"
#include "PETScVectorField.h"

//#include "PETScForceDipole.h"
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
typedef GlobalProperty<1,2> GlobalVector2D;   
typedef GlobalProperty<1,3> GlobalVector3D;   
typedef GlobalProperty<4,2> NonAffineTensor2D;   
typedef GlobalProperty<4,3> NonAffineTensor3D;   

PYBIND11_MODULE(_elasticity, m)
{
    int external_init = slepc::initialize();

    if (!external_init)
    {
        Py_AtExit(slepc::finalize);
    }
    //Export PETSc Managers
    export_EigenManager(m);
    export_PETScManager(m);
    
    //Export PETScObservables
    export_PETScVectorFieldBase(m);
    export_PETScVectorField< PETScVectorField2D >(m,"PETScVectorField2D");
    export_PETScVectorField< PETScVectorField3D >(m,"PETScVectorField3D");
    
    export_GlobalPropertyBase(m);
    export_GlobalProperty< GlobalVector2D >(m,"GlobalVector2D");    
    export_GlobalProperty< GlobalVector3D >(m,"GlobalVector3D");    
    export_GlobalProperty< NonAffineTensor2D >(m,"NonAffineTensor2D");    
    export_GlobalProperty< NonAffineTensor3D >(m,"NonAffineTensor3D");    
    export_GlobalScalar(m,"GlobalScalar");    
    
    //Export PETSc Hessian Objects
    export_HessianBase(m);
    export_PETScHessianBase(m);
    //export_EigenHessianBase(m);
    //export_PETScHessian< PETScHessian<2> >(m, "PETScHessian2D");
    //export_PETScHessian< PETScHessian<3> >(m, "PETScHessian3D");
    export_SLEPcHessian< SLEPcHessian<2> >(m, "SLEPcHessian2D");
    export_SLEPcHessian< SLEPcHessian<3> >(m, "SLEPcHessian3D");
    //export_SpectraHessian< SpectraHessian<2> >(m, "SpectraHessian2D");
    
    //Export PETSc Calculators
    //export_PETScForceDipoleCalculator(m);
    export_SLEPcNMA(m);
    //export_SLEPcHessian(m);
    //export_SpectraHessian(m);
}
