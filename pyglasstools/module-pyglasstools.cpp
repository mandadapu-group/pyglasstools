#include "SimBox.h"
#include "ParticleSystem.h"
#include "Calculator.h"
#include "Manager.h"
#include "MPIInterface.h"
//#include "MPIBasicOperations.h"

#include "extern/pybind11/include/pybind11/pybind11.h"

//Typedefs several pair potentials here

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the pyglasstools python module and define the exports here.
 */


namespace MPI
{
    //! Initialize the MPI environment
    int initialize()
    {
        // initialize MPI if it has not been initialized by another program
        int external_init = 0;
        MPI_Initialized(&external_init);
        if (!external_init)
            {
                MPI_Init(0, (char ***) NULL);
            }

        return external_init;
    }

    //! Get the processor name associated to this rank
    std::string get_processor_name()
    {
        char proc_name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(proc_name, &name_len);
        return std::string(proc_name);
    }

    //! Finalize MPI environment
    void finalize()
    {
        MPI_Finalize();
    }
}

PYBIND11_MODULE(_pyglasstools, m)
{
    int external_init = MPI::initialize();

    // if HOOMD called MPI_Init, it should call MPI_Finalize at exit
    if (!external_init)
    {
        Py_AtExit(MPI::finalize);
    }
    export_SimBox(m);
    
    export_ParticleSystem(m);
    
    export_GlobalCalculator(m);    
    
    export_LocalCalculator(m); 
    
    export_Manager(m); 
    
    export_MPICommunicator(m);  
}
