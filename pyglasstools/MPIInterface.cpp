// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file was part of the HOOMD-blue project, released under the BSD 3-Clause License.

#include "MPIInterface.h"

#include <stdexcept>
#include <iostream>
#include <sstream>

/*! \file MPI::Communicator.cc
    \brief Defines MPI::Communicator and related classes
*/

//! Default constructor
MPI::Communicator::Communicator()
    : m_rank(0), m_n_rank(1)
{
    m_comm = m_comm_world = MPI_COMM_WORLD;

    // use all ranks in a single partition
    int size;
    MPI_Comm_size(m_comm_world, &size);
    m_n_rank = size;

    int rank;
    MPI_Comm_rank(m_comm, &rank);
    m_rank = rank;
}

//! Returns the MPI communicator
MPI_Comm MPI::Communicator::getCommunicator() const
{
    return m_comm;
}

//! Returns the World MPI communicator
MPI_Comm MPI::Communicator::getWorldCommunicator() const
{
    return m_comm_world;
}

//!< Partition the communicator
/*! \param nrank Number of ranks per partition
*/
void MPI::Communicator::splitPartitions(unsigned int nrank)
{
    int num_total_ranks;
    MPI_Comm_size(m_comm_world, &num_total_ranks);

    unsigned int partition = 0;
    m_n_rank = nrank;

    if (m_n_rank == 0)
        throw std::runtime_error("--nrank setting has to be > 0");

    int rank;
    MPI_Comm_rank(m_comm_world, &rank);

    if (num_total_ranks % m_n_rank != 0)
        throw std::runtime_error("Invalid setting --nrank.");

    partition = rank / m_n_rank;

    // Split the communicator
    MPI_Comm new_comm;
    MPI_Comm_split(m_comm_world, partition, rank, &new_comm);

    // update communicator
    m_comm = new_comm;

    MPI_Comm_rank(m_comm, &rank);
    m_rank = rank;
}

//! Return the rank of this processor in the partition
unsigned int MPI::Communicator::getRank() const
{
    return m_rank;
}

//! Return the global rank of this processor
unsigned int MPI::Communicator::getRankGlobal() const
{
    // get rank on world communicator
    int rank;
    MPI_Comm_rank(m_comm_world, &rank);
    return rank;
}

//! Return the number of ranks in this partition
unsigned int MPI::Communicator::getNRanks() const
{
    int size;
    MPI_Comm_size(m_comm, &size);
    return size;
}

//! Return the global communicator size
unsigned int MPI::Communicator::getSizeGlobal() const
{
    int size;
    MPI_Comm_size(m_comm_world, &size);
    return size;
}

//! Returns the partition number of this processor
unsigned int MPI::Communicator::getPartition() const
{
    return getRankGlobal()/m_n_rank;
}

//! Returns the number of partitions
unsigned int MPI::Communicator::getNPartitions() const
{
    return getSizeGlobal()/m_n_rank;
}

//! Returns true if this is the root processor
bool MPI::Communicator::isRoot() const
{
    return getRank() == 0;
}

//! Perform a job-wide MPI barrier
void MPI::Communicator::barrier()
{
    MPI_Barrier(m_comm);
}

//! Perform a job-wide MPI Abort
void MPI::Communicator::mpiabort(int error_code)
{
    if( getSizeGlobal() > 1)
    {
        // delay for a moment to give time for error messages to print
        int msec = 1000;
        usleep(msec*1000);
        MPI_Abort(m_comm, error_code);
    }
}

void export_MPICommunicator(pybind11::module& m)
{
    pybind11::class_<MPI::Communicator, std::shared_ptr<MPI::Communicator> > mpicommunicator(m,"Communicator");
    mpicommunicator.def(pybind11::init< >())
    .def("splitPartitions", &MPI::Communicator::splitPartitions)
    .def("getPartition", &MPI::Communicator::getPartition)
    .def("getNRanks", &MPI::Communicator::getNRanks)
    .def("getRank", &MPI::Communicator::getRank)
    .def("barrier", &MPI::Communicator::barrier)
    .def("abort", &MPI::Communicator::mpiabort)
    .def("scatter_v", &MPI::Communicator::pyscatter_v< std::vector< Eigen::MatrixXd > >)
    .def("all_gather_v", &MPI::Communicator::pyall_gather_v< std::vector< Eigen::MatrixXd > >)
    .def("scatter_v", &MPI::Communicator::pyscatter_v< std::vector< double > >)
    .def("all_gather_v", &MPI::Communicator::pyall_gather_v< std::vector< double > >)
    .def("bcast", &MPI::Communicator::pybcast< double >)
    .def("bcast", &MPI::Communicator::pybcast< int >)
    .def("send", &MPI::Communicator::send< bool >)
    .def("recv", &MPI::Communicator::pyrecv< bool >)
    .def("getSizeGlobal", &MPI::Communicator::getSizeGlobal)
    .def("getRankGlobal", &MPI::Communicator::getRankGlobal)
    ;
}
