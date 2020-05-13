// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

// ensure that HOOMDMath.h is the first thing included
#include <mpi.h>

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
    class PYBIND11_EXPORT Communicator
        {
        public:
            //! Constructor with externally provided MPI communicator (only in MPI enabled builds)
            /*! \param mpi_comm world MPI communicator
             */
            Communicator();
            
            //! Destructor
            virtual ~Communicator() {};

            //! Returns the MPI communicator
            MPI_Comm getCommunicator() const
            {
                return m_comm;
            }
            //! Returns the World MPI communicator
            //MPI_Comm getHOOMDWorldCommunicator() const
            //    {
            //    return m_comm_world;
            //    }

            //!< Partition the communicator
            /*! \param nrank Number of ranks per partition
            */
            void splitPartitions(unsigned int nrank);

            //! Return the rank of this processor in the partition
            unsigned int getRank() const
            {
                return m_rank;
            }

            //! Return the global rank of this processor
            unsigned int getRankGlobal() const
            {
                // get rank on world communicator
                int rank;
                MPI_Comm_rank(m_comm_world, &rank);
                return rank;
            }

            //! Return the global communicator size
            unsigned int getNRanksGlobal() const
            {
                int size;
                MPI_Comm_size(m_comm_world, &size);
                return size;
            }

            //! Returns the partition number of this processor
            unsigned int getPartition() const
            {
                return getRankGlobal()/m_n_rank;
            }

            //! Returns the number of partitions
            unsigned int getNPartitions() const
            {
                return getNRanksGlobal()/m_n_rank;
            }

            //! Return the number of ranks in this partition
            unsigned int getNRanks() const;

            //! Returns true if this is the root processor
            bool isRoot() const
            {
                return getRank() == 0;
            }

            //! Perform a job-wide MPI barrier
            void barrier()
            {
                MPI_Barrier(m_comm);
            }

        protected:
            MPI_Comm m_comm;                   //!< The MPI communicator
            MPI_Comm m_comm_world;                  //!< The world communicator
            
            unsigned int m_rank;                   //!< Rank of this processor (0 if running in single-processor mode)
            unsigned int m_n_rank;                 //!< Ranks per partition
        };
}

//! Exports Communicator to python
void export_MPICommunicator(pybind11::module& m);
