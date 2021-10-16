#include "GlobalProperty.h"
#include "HessianBase.h"

#include <pybind11/pybind11.h>

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the observables python module and define the exports here.
 */

typedef GlobalProperty<1,2> GlobalVector2D;   
typedef GlobalProperty<1,3> GlobalVector3D;   
typedef GlobalProperty<4,2> NonAffineTensor2D;   
typedef GlobalProperty<4,3> NonAffineTensor3D;   

PYBIND11_MODULE(_elasticity, m)
{
    export_GlobalPropertyBase(m);
    export_GlobalProperty< GlobalVector2D >(m,"GlobalVector2D");    
    export_GlobalProperty< GlobalVector3D >(m,"GlobalVector3D");    
    export_GlobalProperty< NonAffineTensor2D >(m,"NonAffineTensor2D");    
    export_GlobalProperty< NonAffineTensor3D >(m,"NonAffineTensor3D");    
    export_GlobalScalar(m,"GlobalScalar");    
    
    export_HessianBase(m);
}
