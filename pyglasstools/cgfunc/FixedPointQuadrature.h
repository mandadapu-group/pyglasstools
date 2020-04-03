//A one-stop wrapper for Quadrature calculations 
//https://scicomp.stackexchange.com/questions/20786/c-library-for-numerical-intergration-quadrature
#ifndef __FIXEDPOINT_QUADRATURE_H__
#define __FIXEDPOINT_QUADRATURE_H__

#include <iostream>
#include <cmath>

#include <functional>
#include <memory>
#include <utility>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

template < typename F >
class gsl_fixedpoint_quad
{
    protected:
        F f;
        int limit;
        std::unique_ptr < gsl_integration_glfixed_table,
                        std::function < void(gsl_integration_glfixed_table*) >
                        > workspace;

        static double gsl_wrapper(double x, void * p)
        {
            gsl_fixedpoint_quad * t = reinterpret_cast<gsl_fixedpoint_quad*>(p);
            return t->f(x);
        }

    public:
        gsl_fixedpoint_quad(F f, int limit)
            : f(f)
            , limit(limit)
            , workspace(gsl_integration_glfixed_table_alloc(limit), gsl_integration_glfixed_table_free)
        {}
        virtual ~gsl_fixedpoint_quad(){};

        double integrate(double min, double max)
        {
            gsl_function gsl_f;
            gsl_f.function = &gsl_wrapper;
            gsl_f.params = this;
            
            return gsl_integration_glfixed ( &gsl_f, min, max, workspace.get());
        }
};

template < typename F >
double GSLFixedPointQuadrature(F func,
            std::pair<double,double> const& range,
            int limit = 5)
{
    return gsl_fixedpoint_quad<F>(func, limit).integrate(range.first, range.second);
}

#endif //__FIXEDPOINT_QUADRATURE_H__
