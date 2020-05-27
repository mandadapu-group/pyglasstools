#ifndef __PETSCLINEARSOLVER_H__
#define __PETSCLINEARSOLVER_H__

#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>
#include <memory>

#include <pyglasstools/MathAndTypes.h>
#include <pyglasstools/ParticleSystem.h>
#include <pyglasstools/potential/PairPotential.h>
#include "NonAffineManager.h"

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
namespace py = pybind11;

#include <Aboria.h>
namespace abr = Aboria;

#include <petscksp.h>


//A Class to solve A*x=b problems
class ForceDipoleSolver
{
    private:
        std::shared_ptr< Mat > m_hessian;
        std::shared_ptr< ParticleSystem > m_sysdata;
        Vec x;
        Vec b;
        KSP ksp;
    public:
        ForceDipoleSolver(std::shared_ptr< ParticleSystem > sysdata, Mat hessian)
            : m_sysdata(sysdata), m_hessian(std::make_shared<Mat>(hessian))
        {
        }
};
#endif
