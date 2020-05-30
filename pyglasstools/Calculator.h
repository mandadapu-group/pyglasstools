#ifndef __CALCULATOR_H__
#define __CALCULATOR_H__

#include <vector>
#include <cmath>
#include <string>
#include <map>

#include <algorithm>
#include <memory>

#include "MathAndTypes.h"
#include "ParticleSystem.h"
#include "potential/PairPotential.h"
#include "Observable.h"

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <Aboria.h>
namespace abr = Aboria;

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Core>

//Public interface for classes computing observables
class PYBIND11_EXPORT Calculator
{
    protected:

        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< PairPotential > m_potential;
        double max_rcut;
    public:
        
        Calculator( std::shared_ptr< ParticleSystem > sysdata, 
                    std::shared_ptr< PairPotential > potential)
            : m_sysdata(sysdata), m_potential(potential)
            {
                double maxforce_rcut = potential->scaled_rcut;
                double maxdiameter = *std::max_element( std::begin(abr::get<diameter>(m_sysdata->particles)), 
                                                        std::end(abr::get<diameter>(m_sysdata->particles)) );
                maxforce_rcut *= maxdiameter;
                max_rcut = std::max(maxforce_rcut,maxdiameter); //whichever is maximum
            };
        ~Calculator(){};
        

        void setSystemData( std::shared_ptr< ParticleSystem > sysdata)
        {
            m_sysdata = sysdata;
        } 
        
        virtual void compute(){};
};

void export_Calculator(py::module& m)
{
    py::class_<Calculator, std::shared_ptr<Calculator> >(m,"Calculator")
    .def(py::init< std::shared_ptr< ParticleSystem >, std::shared_ptr< PairPotential >  >())
    ;
};

#endif //__CALCULATOR_H__
