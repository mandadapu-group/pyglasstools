#include "SimBox.h"
#include "ParticleSystem.h"
#include "Calculator.h"
#include "Manager.h"
#include "MPIInterface.h"
#include "MPIBasicOperations.h"

#include "extern/pybind11/include/pybind11/pybind11.h"

//Typedefs several pair potentials here

//! Create the python module
/*! each class setup their own python exports in a function export_ClassName
 create the pyglasstools python module and define the exports here.
 */

//void mpi_barrier_world()
//  {
//   #ifdef ENABLE_MPI
//  MPI_Barrier(MPI_COMM_WORLD);
//  #endif
//  }

// values used in measuring hoomd launch timing
unsigned int launch_time, start_time, mpi_init_time;
bool launch_timing=false;



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
    
    //export_MPICommunicator(m);  
}
