#include "CoarseGrainFunction.h"
#include "CGFuncOctic.h"
#include "CGFuncRect.h"
#include "CGFuncMollifier.h"

#include "../extern/pybind11/include/pybind11/pybind11.h"
namespace py = pybind11;

//Typedefs several coarsegrain functions  here
typedef PolynomialCGFunc<Octic> CGFuncOctic;
typedef PolynomialCGFunc<Rect> CGFuncRect;
typedef GeneralCGFunc<Mollifier> CGFuncMollifier;

PYBIND11_MODULE(_cgfunc, m)
{
    export_CoarseGrainFunction(m);
    export_PolynomialCGFunc<CGFuncOctic>(m, "CGFuncOctic");
    export_PolynomialCGFunc<CGFuncRect>(m, "CGFuncRect");
    export_GeneralCGFunc<CGFuncMollifier>(m, "CGFuncMollifier");
}
