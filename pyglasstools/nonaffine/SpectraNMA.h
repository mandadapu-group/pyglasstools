#ifndef __SPECTRA_NMA_H__
#define __SPECTRA_NMA_H__

#include "HessianBase.h"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

class PYBIND11_EXPORT SpectraNMA
{
    protected:
        std::shared_ptr< PETScHessianBase > m_hessian; 
        
        //List of observables to be computed!
        std::map<std::string, std::shared_ptr< PETScVectorFieldBase > > m_vectorfields;
        std::map<std::string, std::shared_ptr< PETScGlobalPropertyBase > > m_observables;
        
        int nconv; 
        double m_maxeigval;
        double m_mineigval;
        
        //PETSc objects for all computation 
        PetscErrorCode ierr; 
    
    public:
        SLEPcNMA(std::shared_ptr< PETScHessianBase > hessian);
        virtual ~SLEPcNMA()
        {
        };
        void setHessian(std::shared_ptr< PETScHessianBase > hessian)
        {
            m_hessian = hessian;
        }
        virtual void addVectorField(const std::shared_ptr< PETScVectorFieldBase >& obs)
        {
            m_vectorfields.insert(std::pair<std::string, std::shared_ptr< PETScVectorFieldBase > >(obs->name, obs));
        }

        virtual void addGlobalProperty(const std::shared_ptr< PETScGlobalPropertyBase >& obs)
        {
            m_observables.insert(std::pair<std::string, std::shared_ptr< PETScGlobalPropertyBase > >(obs->name, obs));
        }
        
        void getAllEigenPairs(std::string package);
        void getEigenPairs(std::string package);
        //Another function for getting some eigenpairs
        void getMaxEigenvalue();

        void saveResults(EPS& eps)
        {
            if (m_observables.count("nconv") > 0)
            {
                m_observables["nconv"]->setValue((double)nconv);
            }
            //Save Solution
            if(nconv >0)
            {
                calculateNonAffineTensor(eps);

                m_hessian->m_manager->printPetscNotice(5,"Storing selected eigenpairs \n");
                
                std::string prefix = "eigenvector_";
                std::string prefix1 = "eigenvalue_";
                std::string prefix2 = "eigenrelerror_";
                std::stringstream filename;
                for(int i = 0; i < nconv; ++i)
                {
                    filename.str("");
                    filename << prefix << i;
                    if (m_vectorfields.count(filename.str()) > 0)
                    {
                        m_vectorfields[filename.str()]->createVector(m_hessian->hessian);
                        ierr = EPSGetEigenvector(eps,i,m_vectorfields[filename.str()]->vectorobs,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                    }

                    filename.str("");
                    filename << prefix1 << i;
                    if (m_observables.count(filename.str()) > 0)
                    {
                        PetscReal vr;
                        ierr = EPSGetEigenvalue(eps,i,&vr,NULL);CHKERRABORT(PETSC_COMM_WORLD,ierr);
                        m_observables[filename.str()]->setValue(vr);
                    }
                    
                    filename.str("");
                    filename << prefix2 << i;
                    if (m_observables.count(filename.str()) > 0)
                    {
                        PetscReal error;
                        ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error); CHKERRABORT(PETSC_COMM_WORLD,ierr);
                        m_observables[filename.str()]->setValue(error);
                    }
                }
            }
        }
        
        //void getEigenPairs();
        void calculateNonAffineTensor(const EPS& eps); 
        
        //double getEigenvalue(unsigned int index); 
        //std::vector<double> getEigenvector(unsigned int index, bool forall); 

        //void saveNonAffineTensor(unsigned int i, unsigned int j, std::shared_ptr<MPI::LogFile> logfile); 
};

SLEPcNMA::SLEPcNMA( std::shared_ptr< PETScHessianBase > hessian) 
    : m_hessian(hessian), nconv(0), m_maxeigval(std::numeric_limits<double>::max()), m_mineigval(-std::numeric_limits<double>::max())
{
    ierr = PetscOptionsInsertString(NULL,m_hessian->m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
};

void SLEPcNMA::getMaxEigenvalue()
{}

void SLEPcNMA::getEigenPairs(std::string package)
{}

void SLEPcNMA::getAllEigenPairs(std::string package)
{}

void SLEPcNMA::calculateNonAffineTensor(const EPS& eps) 
{}

void export_SLEPcNMA(py::module& m)
{
    py::class_<SLEPcNMA, std::shared_ptr<SLEPcNMA> >(m,"SLEPcNMA")
    .def(py::init< std::shared_ptr< PETScHessianBase > >())
    .def("setHessian", &SLEPcNMA::setHessian)
    .def("addVectorField", &SLEPcNMA::addVectorField)
    .def("addGlobalProperty", &SLEPcNMA::addGlobalProperty)
    .def("getAllEigenPairs", &SLEPcNMA::getAllEigenPairs)
    //.def("calculateNonAffineTensor", &T::calculateNonAffineTensor)
    //.def("getEigenPairs", &T::getEigenPairs)
    //.def("saveNonAffineTensor", &T::saveNonAffineTensor)
    //.def("saveEigenvector", &T::saveEigenvector)
    //.def("getEigenvalue", &T::getEigenvalue)
    //.def("saveEigenvalue", &T::saveEigenvalue)
    //.def("getPETScMatRange", &T::getPETScMatRange)
    ;
};

#endif //__SPECTRA_NMA_H__
