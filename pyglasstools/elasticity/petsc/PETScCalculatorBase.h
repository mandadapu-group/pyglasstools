#ifndef __PETSC_CALCULATOR_BASE_H__
#define __PETSC_CALCULATOR_BASE_H__

#include "PETScHessianBase.h"

/*
 * A base class for constructing calculators that utilize PETSc. This may include eigensolvers or linearsolvers. 
 * The goal is to provide common interface and data between different calculators such as the ability to save and store observables as well as information about null space. 
 */
class PYBIND11_EXPORT PETScCalculatorBase
{
    protected:

        std::shared_ptr< PETScHessianBase > m_hessian; //<! pointer to the Hessian base class.
        
        std::map<std::string, std::shared_ptr< PETScVectorFieldBase > > m_vectorfields; //<! list of vector-field observables to compute
        std::map<std::string, std::shared_ptr< GlobalPropertyBase > > m_observables; //<! list of global observables to compute
        //PETSc error code 
        PetscErrorCode ierr; 
        
        MatNullSpace nullspace; //!< a class for storing the null space of the Hessian, i.e., spanned by the zero modes.
        Vec* transmodes; 
    public:
        PETScCalculatorBase(std::shared_ptr< PETScHessianBase > hessian);
        virtual ~PETScCalculatorBase()
        {
            int Dim = m_hessian->m_sysdata->simbox->dim;
            for(int i = 0; i < Dim; ++i)
            {
                VecDestroy(&transmodes[i]);
            }
            PetscFree(transmodes);
            MatNullSpaceDestroy(&nullspace);
        };
        
        /* Set the stored Hessian with another Hessian. Args:
         * hessian: a Hessian class 
         */
        void setHessian(std::shared_ptr< PETScHessianBase > hessian)
        {
            m_hessian = hessian;
            MatSetNullSpace(m_hessian->hessian,nullspace);
        }

        /* Add a new vector field observable to a list of existing ones. Args:
         * obs: the new vector field observable
         */
        virtual void addVectorField(const std::shared_ptr< PETScVectorFieldBase >& obs)
        {
            m_vectorfields.insert(std::pair<std::string, std::shared_ptr< PETScVectorFieldBase > >(obs->name, obs));
        }
        
        /* Add a new global observable/property toa  list of existing ones. Args:
         * obs: the new global observable
         */
        virtual void addGlobalProperty(const std::shared_ptr< GlobalPropertyBase >& obs)
        {
            m_observables.insert(std::pair<std::string, std::shared_ptr< GlobalPropertyBase > >(obs->name, obs));
        }

        void setNullSpaceBasis(PetscInt id_i, PetscInt Istart, PetscInt Iend, PetscInt real_id);
};

/* Parametrized constructor. Args:
 * hessian: the input Hessian matrix
 */
PETScCalculatorBase::PETScCalculatorBase( std::shared_ptr< PETScHessianBase > hessian) 
    : m_hessian(hessian)
{
    int Dim = m_hessian->m_sysdata->simbox->dim;
    //Add any relevant command-line options for the PETSc linear algebra solver.
    ierr = PetscOptionsInsertString(NULL,m_hessian->m_manager->cmd_line_options.c_str());CHKERRABORT(PETSC_COMM_WORLD,ierr);
    PetscInt Istart; PetscInt Iend; MatGetOwnershipRange(m_hessian->hessian, &Istart, &Iend);
    m_hessian->m_manager->printPetscNotice(5,"Begin Assembling the null space of PETSc matrix\n");
    PetscMalloc1(Dim,&transmodes);    
    for (int i = 0; i < Dim; ++i)
    {
        MatCreateVecs(m_hessian->hessian,NULL,&transmodes[i]);
    }
    
    for (PetscInt i = Istart; i < Iend; ++i) 
    {
        //Instantiate the iterator at a particular position
        auto p_i = m_hessian->m_sysdata->particles.begin()+(int)(i/2);
        PetscInt id_i = abr::get<abr::id>(*p_i);
        setNullSpaceBasis(id_i, Istart, Iend, i);
    }
    
    m_hessian->m_manager->printPetscNotice(5,"Assemble the null space of PETSc matrix\n");
    for (int i = 0; i < Dim; ++i)
    {
        VecAssemblyBegin(transmodes[i]);
        VecAssemblyEnd(transmodes[i]);
        VecNormalize(transmodes[i],NULL);
    }
    MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_FALSE,Dim,transmodes,&nullspace);
    MatSetNullSpace(m_hessian->hessian,nullspace);
};

void PETScCalculatorBase::setNullSpaceBasis(PetscInt id_i, PetscInt Istart, PetscInt Iend, PetscInt real_id)
{
    int Dim = m_hessian->m_sysdata->simbox->dim;
    if (Dim*id_i == real_id)
    { 
        VecSetValue(transmodes[0],Dim*id_i,1, INSERT_VALUES);
        VecSetValue(transmodes[1],Dim*id_i,0, INSERT_VALUES);
        if (Dim == 3)
            VecSetValue(transmodes[Dim-1],Dim*id_i,0, INSERT_VALUES);
    }
    //y-component of the row
    else if (Dim*id_i+1 == real_id)
    {
        VecSetValue(transmodes[0],Dim*id_i+1,0, INSERT_VALUES);
        VecSetValue(transmodes[1],Dim*id_i+1,1, INSERT_VALUES);
        if (Dim == 3)
            VecSetValue(transmodes[Dim-1],Dim*id_i+1,0, INSERT_VALUES);
    }
    
    //z-component of the row
    else if (Dim*id_i+2 == real_id && Dim == 3)
    {
        VecSetValue(transmodes[0],Dim*id_i+2,0, INSERT_VALUES);
        VecSetValue(transmodes[1],Dim*id_i+2,0, INSERT_VALUES);
        if (Dim == 3)
            VecSetValue(transmodes[Dim-1],Dim*id_i+2,1, INSERT_VALUES);
    }
};

/*
 * Helper function to exprot PETScCalculatorBase to Python
 */
void export_PETScCalculatorBase(py::module& m)
{
    py::class_<PETScCalculatorBase, std::shared_ptr<PETScCalculatorBase> >(m, "PETScCalculatorBase")
    .def(py::init< std::shared_ptr< PETScHessianBase > >())
    .def("setHessian", &PETScCalculatorBase::setHessian)
    .def("addVectorField", &PETScCalculatorBase::addVectorField)
    .def("addGlobalProperty", &PETScCalculatorBase::addGlobalProperty)
    ;
};

#endif //__PETSC_CALCULATOR_BASE_H__
