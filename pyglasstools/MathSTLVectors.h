#ifndef __VECTOR_MATH_H__
#define __VECTOR_MATH_H__

#include <algorithm>
#include <functional>
#include <vector>
#include <cassert>

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

struct SymmetricTensorField
{
    std::vector<double> XX;
    std::vector<double> YY;
    std::vector<double> XY;
    std::vector<double> ZZ;
    std::vector<double> XZ;
    std::vector<double> YZ;

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

typedef std::vector<double> ScalarField;

#endif
