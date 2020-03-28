#include "Quadrature.h"
#include "CoarseGrainFunction.h"
#include "CGFuncOctic.h"

#include "../extern/pybind11/include/pybind11/pybind11.h"
namespace py = pybind11;

//Typedefs several coarsegrain functions  here
typedef ShortRangeCGFunc<Octic> CGFuncOctic;

PYBIND11_MODULE(_cgfunc, m)
{
    export_CoarseGrainFunction(m);
    export_ShortRangeCGFunc<CGFuncOctic>(m, "CGFuncOctic");
}
