#ifndef __VECTOR_MATH_H__
#define __VECTOR_MATH_H__

#include <algorithm>
#include <functional>
#include <vector>
#include <cassert>
#include <Eigen/Dense>

#include <Aboria.h>

ABORIA_VARIABLE(velocity, Aboria::vdouble3, "velocity");
ABORIA_VARIABLE(displacement, Aboria::vdouble3, "displacement");
ABORIA_VARIABLE(mass, double, "mass");
ABORIA_VARIABLE(diameter, double, "diameter");
typedef Aboria::Particles< std::tuple<velocity, displacement,diameter, mass> > AboriaParticles;
typedef typename AboriaParticles::position position;        

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
};

template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::minus<T>());
    return result;
};

template <typename T>
std::vector<T> operator*(const std::vector<T>& a, const T& b)
{
    std::vector<T> result;
    result = a;//.reserve(a.size());
    std::transform( result.begin(), result.end(), result.begin(),
                    std::bind(  std::multiplies<T>(), 
                                std::placeholders::_1, b)
                    );
    return result;
};

template <typename T>
std::vector<T> operator*(const T& b, const std::vector<T>& a)
{
    std::vector<T> result;
    result = a;//.reserve(a.size());
    std::transform( result.begin(), result.end(), result.begin(),
                    std::bind(  std::multiplies<T>(), 
                                std::placeholders::_1, b)
                    );
    return result;
};

struct TensorField
{
    std::vector<double> XX;
    std::vector<double> YY;
    std::vector<double> XY;
    std::vector<double> YX;
    std::vector<double> ZZ;
    std::vector<double> XZ;
    std::vector<double> ZX;
    std::vector<double> YZ;
    std::vector<double> ZY;

    void resize(unsigned int size)
    {
        XX.resize(size); 
        YY.resize(size); 
        XY.resize(size); 
        ZZ.resize(size); 
        XZ.resize(size); 
        XX.resize(size); 
    }
};


typedef Eigen::MatrixXd Tensor;//ScalarField;

struct VectorField
{
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;

    void resize(unsigned int size)
    {
        X.resize(size); 
        Y.resize(size); 
        Z.resize(size); 
    }
};

struct Vector
{
    double X;
    double Y;
    double Z;
};

typedef std::vector<double> ScalarField;
typedef double Scalar;

#endif
