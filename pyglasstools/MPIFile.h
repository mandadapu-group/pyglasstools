#ifndef __MPI_FILE_H__
#define __MPI_FILE_H__

// ensure that HOOMDMath.h is the first thing included
#include <sstream>
#include <vector>
#include "MPIInterface.h"
/*! \file Communicator.h
    \brief Declares Communicator, which initializes the MPI environment
*/
#include "extern/pybind11/include/pybind11/pybind11.h"

//! Defines the MPI configuration for the simulation
/*! \ingroup data_structs
    Communicator is class that stores the MPI communicator and splits it into partitions if needed.

    It is constructed *before* ExecutionConfiguration and other classes (Messenger) that need the MPI
    world communicator.
*/
namespace MPI
{

    class PYBIND11_EXPORT LogFile
    {
        private:
            MPI_File m_file;
            std::shared_ptr<MPI::Communicator> m_comm;
        public:
            //! Constructor with externally provided MPI communicator (only in MPI enabled builds)
            /*! \param mpi_comm world MPI communicator
             */
            LogFile(std::shared_ptr<MPI::Communicator> comm, std::string filename)
                : m_comm(comm)
            {
                MPI_File_open(m_comm->getCommunicator(), filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &m_file); 
                MPI_File_set_atomicity(m_file, true);
            };
            
            //! Destructor
            virtual ~LogFile()
            {
                MPI_File_sync(m_file);
                MPI_File_close(&m_file);
            };

            void write_shared(std::string input)
            {
                //serialize input?
                MPI_File_write_shared(m_file, input.c_str(), input.size(), MPI_CHAR, NULL);
            };
            void write_ordered(std::string input)
            {
                //serialize input?
                MPI_File_write_ordered(m_file, input.c_str(), input.size(), MPI_CHAR, NULL);
            };
    };
}

//! Exports Communicator to python
void export_MPIFile(pybind11::module& m)
{
    pybind11::class_<MPI::LogFile, std::shared_ptr<MPI::LogFile> > mpifile(m,"MPILogFile");
    mpifile.def(pybind11::init< std::shared_ptr<MPI::Communicator>, std::string >())
    .def("write_shared", &MPI::LogFile::write_shared)
    ;
}

#endif
