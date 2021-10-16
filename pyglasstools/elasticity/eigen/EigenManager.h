#ifndef __EIGEN_MANAGER_H__
#define __EIGEN_MANAGER_H__

#include <pyglasstools/Manager.h>
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <pybind11/embed.h>
namespace py = pybind11;


class PYBIND11_EXPORT EigenManager : public Manager
{
    public:
        double pinv_tol, eigen_tol;  
        int nev, ncv, maxiter;
        std::string selrule;     
        double fd_random_min, fd_random_max; 
        std::string fd_mode; 
        
        EigenManager() : pinv_tol(1e-12), eigen_tol(1e-6), nev(1), ncv(2), maxiter(100), selrule("LM"),fd_random_min(0), fd_random_max(std::numeric_limits<double>::max()), fd_mode("random")
        {
            detail::argparser<double>("-pinv_tol",pinv_tol, "[ERROR] Invalid value for pinv_tol", argv);
            detail::argparser<double>("-eigen_tol",eigen_tol, "[ERROR] Invalid value for lowerbound_tol", argv);
            detail::argparser<int>("-spectra_nev",nev, "[ERROR] Invalid value spectra_nev", argv);
            detail::argparser<int>("-spectra_ncv",ncv, "[ERROR] Invalid value spectra_ncv", argv);
            detail::argparser<int>("-spectra_maxiter",maxiter, "[ERROR] Invalid value for maxier", argv);
            detail::argparser<std::string>("-spectra_selrule",selrule, "[ERROR] Invalid value for spectra selrule", argv);
        };
        ~EigenManager(){};
};


void export_EigenManager(py::module& m)
{
    py::class_<EigenManager, std::shared_ptr<EigenManager> >(m,"EigenManager")
    .def(py::init<>())
    ;
};

#endif
