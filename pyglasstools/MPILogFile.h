#ifndef __MPI_LOGFILE_H__
#define __MPI_LOGFILE_H__

// ensure that HOOMDMath.h is the first thing included
#include <sstream>
#include <vector>
#include "MPICommunicator.h"
#include "LogFile.h"
#include <pybind11/pybind11.h>

namespace MPI
{

    class PYBIND11_EXPORT ParallelLogFile : public BaseLogFile
    {
        private:
            MPI_File m_file;
            std::shared_ptr<MPI::ParallelCommunicator> m_comm;
        
        public:
            //! Constructor with externally provided MPI communicator (only in MPI enabled builds)
            /*! \param mpi_comm world MPI communicator
             */
            ParallelLogFile(std::shared_ptr<MPI::ParallelCommunicator> comm, std::string filename)
                : m_comm(comm)
            {
                MPI_File_open(m_comm->getCommunicator(), filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &m_file); 
                MPI_File_set_atomicity(m_file, true);
            };
            
            //! Destructor
            virtual ~ParallelLogFile()
            {
                MPI_File_sync(m_file);
                MPI_File_close(&m_file);
            };

            void write(std::string input)
            {
                //serialize input?
                MPI_File_write_shared(m_file, input.c_str(), input.size(), MPI_CHAR, NULL);
            };
            /*
            void write_ordered(std::string input)
            {
                //serialize input?
                MPI_File_write_ordered(m_file, input.c_str(), input.size(), MPI_CHAR, NULL);
            };
            */
    };
}

//! Exports Communicator to python
void export_ParallelLogFile(pybind11::module& m)
{
    pybind11::class_<MPI::ParallelLogFile, BaseLogFile, std::shared_ptr<MPI::ParallelLogFile> > (m,"ParallelLogFile")
    .def(pybind11::init< std::shared_ptr<MPI::ParallelCommunicator>, std::string >())
    .def("write", &MPI::ParallelLogFile::write)
    ;
}

#endif
