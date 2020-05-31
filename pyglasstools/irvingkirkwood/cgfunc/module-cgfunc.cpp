#include "CoarseGrainFunction.h"
#include "CGFuncOctic.h"
#include "CGFuncRect.h"
#include "CGFuncMollifier.h"

#include <pybind11/pybind11.h>
namespace py = pybind11;

//Typedefs several coarsegrain functions  here
typedef FixedPointCGFunc<Octic> CGFuncOcticFixedPoint;
typedef AdaptiveCGFunc<Octic> CGFuncOcticAdaptive;

typedef FixedPointCGFunc<Rect> CGFuncRectFixedPoint;
typedef AdaptiveCGFunc<Rect> CGFuncRectAdaptive;

typedef FixedPointCGFunc<Mollifier> CGFuncMollifierFixedPoint;
typedef AdaptiveCGFunc<Mollifier> CGFuncMollifierAdaptive;

PYBIND11_MODULE(_cgfunc, m)
{
    export_CoarseGrainFunction(m);
    
    export_FixedPointCGFunc<CGFuncOcticFixedPoint>(m, "CGFuncOcticFixedPoint");
    export_AdaptiveCGFunc<CGFuncOcticAdaptive>(m, "CGFuncOcticAdaptive");

    export_FixedPointCGFunc<CGFuncRectFixedPoint>(m, "CGFuncRectFixedPoint");
    export_AdaptiveCGFunc<CGFuncRectAdaptive>(m, "CGFuncRectAdaptive");
    
    export_FixedPointCGFunc<CGFuncMollifierFixedPoint>(m, "CGFuncMollifierFixedPoint");
    export_AdaptiveCGFunc<CGFuncMollifierAdaptive>(m, "CGFuncMollifierAdaptive");
}
