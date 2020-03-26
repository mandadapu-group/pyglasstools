#ifndef __VECTOR_MATH_H__
#define __VECTOR_MATH_H__

#include <algorithm>
#include <functional>
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

#endif
