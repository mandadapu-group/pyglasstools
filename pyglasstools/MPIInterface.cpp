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

unsigned int MPI::Communicator::getNRanks() const
{
    int size;
    MPI_Comm_size(m_comm, &size);
    return size;
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
    .def("getSizeGlobal", &MPI::Communicator::getSizeGlobal)
    .def("getRankGlobal", &MPI::Communicator::getRankGlobal)
    ;
}
