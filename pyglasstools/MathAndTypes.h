#ifndef __VECTOR_MATH_H__
#define __VECTOR_MATH_H__

#include <algorithm>
#include <functional>
#include <vector>
#include <cassert>

#include <Aboria.h>

ABORIA_VARIABLE(velocity, Aboria::vdouble3, "velocity");
ABORIA_VARIABLE(displacement, Aboria::vdouble3, "displacement");
ABORIA_VARIABLE(mass, double, "mass");
ABORIA_VARIABLE(diameter, double, "diameter");
typedef Aboria::Particles< std::tuple<velocity, displacement,diameter, mass> > AboriaParticles;
typedef typename AboriaParticles::position position;        

#endif
