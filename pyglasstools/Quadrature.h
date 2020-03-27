//A one-stop wrapper for Quadrature calculations 
//https://scicomp.stackexchange.com/questions/20786/c-library-for-numerical-intergration-quadrature
#ifndef __QUADRATURE_H__
#define __QUADRATURE_H__

#include <iostream>
#include <cmath>

#include <functional>
#include <memory>
#include <utility>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

template < typename F >
class gsl_quad
{
    protected:
        F f;
        int limit;
        std::unique_ptr < gsl_integration_workspace,
                        std::function < void(gsl_integration_workspace*) >
                        > workspace;

        static double gsl_wrapper(double x, void * p)
        {
            gsl_quad * t = reinterpret_cast<gsl_quad*>(p);
            return t->f(x);
        }

    public:
        gsl_quad(F f, int limit)
            : f(f)
            , limit(limit)
            , workspace(gsl_integration_workspace_alloc(limit), gsl_integration_workspace_free)
        {}
        virtual ~gsl_quad(){};

        double integrate(double min, double max, double epsabs, double epsrel)
        {
            gsl_function gsl_f;
            gsl_f.function = &gsl_wrapper;
            gsl_f.params = this;

            double result, error;
            
            if ( !std::isinf(min) && !std::isinf(max) )
            {
                gsl_integration_qag ( &gsl_f, min, max,
                                     epsabs, epsrel, limit,
                                     6,workspace.get(), &result, &error );
            }
            else if ( std::isinf(min) && !std::isinf(max) )
            {
                gsl_integration_qagil( &gsl_f, max,
                                     epsabs, epsrel, limit,
                                     workspace.get(), &result, &error );
            }
            else if ( !std::isinf(min) && std::isinf(max) )
            {
                gsl_integration_qagiu( &gsl_f, min,
                                     epsabs, epsrel, limit,
                                     workspace.get(), &result, &error );
            }
            else
            {
                gsl_integration_qagi ( &gsl_f,
                                     epsabs, epsrel, limit,
                                     workspace.get(), &result, &error );
            }

            return result;
        }
};

template < typename F >
double GSLQuadrature(F func,
            std::pair<double,double> const& range,
            double epsabs = 1.49e-9, double epsrel = 1.49e-9,
            int limit = 200)
{
    return gsl_quad<F>(func, limit).integrate(range.first, range.second, epsabs, epsrel);
}

#endif //__QUADRATURE_H__
