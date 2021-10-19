#include "MPICommunicator.h"

#include <stdexcept>
#include <iostream>
#include <sstream>

/*! \file MPI::ParallelCommunicator.cc
    \brief Defines MPI::ParallelCommunicator and related classes
*/

//! Default constructor
MPI::ParallelCommunicator::ParallelCommunicator()
    : m_rank(0), m_n_rank(1)
{
    m_comm = MPI_COMM_WORLD;

    // use all ranks in a single partition
    int size;
    MPI_Comm_size(m_comm, &size);
    m_n_rank = size;

    int rank;
    MPI_Comm_rank(m_comm, &rank);
    m_rank = rank;
}

//! Returns the MPI communicator
MPI_Comm MPI::ParallelCommunicator::getCommunicator() const
{
    return m_comm;
}

//! Return the rank of this processor in the partition
unsigned int MPI::ParallelCommunicator::getRank() const
{
    return m_rank;
}

//! Return the number of ranks in this partition
unsigned int MPI::ParallelCommunicator::getNRanks() const
{
    return m_n_rank;
}

//! Returns true if this is the root processor
bool MPI::ParallelCommunicator::isRoot() const
{
    return m_rank == 0;
}

//! Perform a job-wide MPI barrier
void MPI::ParallelCommunicator::barrier()
{
    MPI_Barrier(m_comm);
}

//! Perform a job-wide MPI Abort
void MPI::ParallelCommunicator::abort(int error_code)
{
    // delay for a moment to give time for error messages to print
    int msec = 1000;
    usleep(msec*1000);
    MPI_Abort(m_comm, error_code);
}

void export_ParallelCommunicator(pybind11::module& m)
{
    pybind11::class_<MPI::ParallelCommunicator, std::shared_ptr<MPI::ParallelCommunicator> > (m,"ParallelCommunicator")
    .def(pybind11::init< >())
    .def("getRank", &MPI::ParallelCommunicator::getRank)
    .def("getNRanks", &MPI::ParallelCommunicator::getNRanks)
    .def("isRoot", &MPI::ParallelCommunicator::isRoot)
    .def("barrier", &MPI::ParallelCommunicator::barrier)
    .def("abort", &MPI::ParallelCommunicator::abort)
    .def("scatter_v", &MPI::ParallelCommunicator::pyscatter_v< std::vector< Eigen::MatrixXd > >)
    .def("all_gather_v", &MPI::ParallelCommunicator::pyall_gather_v< std::vector< Eigen::MatrixXd > >)
    .def("scatter_v", &MPI::ParallelCommunicator::pyscatter_v< std::vector< double > >)
    .def("all_gather_v", &MPI::ParallelCommunicator::pyall_gather_v< std::vector< double > >)
    .def("bcast", &MPI::ParallelCommunicator::pybcast< double >)
    .def("bcast", &MPI::ParallelCommunicator::pybcast< int >)
    .def("send", &MPI::ParallelCommunicator::send< bool >)
    .def("recv", &MPI::ParallelCommunicator::pyrecv< bool >)
    ;
}
