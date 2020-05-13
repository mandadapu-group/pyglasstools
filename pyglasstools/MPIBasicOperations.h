// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file was part of the HOOMD-blue project, released under the BSD 3-Clause License.

#ifndef __MPI_INTERFACE_H__
#define __MPI_INTERFACE_H__

/*! \file HOOMDMPI.h
    \brief Defines common functions for MPI operations

    The functions provided here imitate some basic boost.MPI functionality.
*/

#include <mpi.h>

#include <sstream>
#include <vector>

#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp> 
#include <cereal/archives/binary.hpp>


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
                py::print("Test");
                MPI_Init(0, (char ***) NULL);
                py::print("Test 2");
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

    //! Abort MPI runs
    void abort(std::shared_ptr<MPI::Communicator> comm, unsigned int milliseconds)
    {
        if(comm->getNRanksGlobal() > 1)
        {
            // delay for a moment to give time for error messages to print
            usleep(1000 * milliseconds);
            MPI_Abort(comm->getCommunicator(), MPI_ERR_OTHER);
        }
    }
    //! Wrapper around MPI_Bcast that handles any serializable object
    template<typename T>
    void bcast(T& val, unsigned int root, const MPI_Comm mpi_comm)
        {
        int rank;
        MPI_Comm_rank(mpi_comm, &rank);

        char *buf = NULL;
        int recv_count;
        if (rank == (int)root)
            {
            std::stringstream s(std::ios_base::out | std::ios_base::binary);
            cereal::BinaryOutputArchive ar(s);

            // serialize object
            ar << val;

            // do not forget to flush stream
            s.flush();

            // copy string to send buffer
            std::string str = s.str();
            recv_count = str.size();
            buf = new char[recv_count];
            str.copy(buf, recv_count);
            }

        MPI_Bcast(&recv_count, 1, MPI_INT, root, mpi_comm);
        if (rank != (int) root)
            buf = new char[recv_count];

        MPI_Bcast(buf, recv_count, MPI_BYTE, root, mpi_comm);

        if (rank != (int)root)
            {
            // de-serialize
            std::stringstream s(std::string(buf, recv_count), std::ios_base::in | std::ios_base::binary);
            cereal::BinaryInputArchive ar(s);

            ar >> val;
            }

        delete[] buf;
        }

    //! Wrapper around MPI_Scatterv that scatters a vector of serializable objects
    template<typename T>
    void scatter_v(const std::vector<T>& in_values, T& out_value, unsigned int root, const MPI_Comm mpi_comm)
        {
        int rank;
        int size;
        MPI_Comm_rank(mpi_comm, &rank);
        MPI_Comm_size(mpi_comm, &size);

        assert(in_values.size() == (unsigned int) size);

        unsigned int recv_count;
        int *send_counts = NULL;
        int *displs = NULL;

        char *sbuf = NULL;
        if (rank == (int)root)
            {
            send_counts = new int[size];
            displs = new int[size];
            // construct a vector of serialized objects
            typename std::vector<T>::const_iterator it;
            std::vector<std::string> str;
            unsigned int len = 0;
            for (it = in_values.begin(); it!= in_values.end(); ++it)
                {
                unsigned int idx = it - in_values.begin();
                std::stringstream s(std::ios_base::out | std::ios_base::binary);
                cereal::BinaryOutputArchive ar(s);

                // serialize object
                ar << *it;
                s.flush();
                str.push_back(s.str());

                displs[idx] = (idx > 0) ? displs[idx-1]+send_counts[idx-1] : 0;
                send_counts[idx] = str[idx].length();
                len += send_counts[idx];
                }

            // pack vector into send buffer
            sbuf = new char[len];
            for (unsigned int idx = 0; idx < in_values.size(); idx++)
                str[idx].copy(sbuf + displs[idx], send_counts[idx]);
            }

        // scatter sizes of serialized vector elements
        MPI_Scatter(send_counts, 1, MPI_INT, &recv_count, 1, MPI_INT, root, mpi_comm);

        // allocate receive buffer
        char *rbuf = new char[recv_count];

        // scatter actual data
        MPI_Scatterv(sbuf, send_counts, displs, MPI_BYTE, rbuf, recv_count, MPI_BYTE, root, mpi_comm);

        // de-serialize data
        std::stringstream s(std::string(rbuf, recv_count), std::ios_base::in | std::ios_base::binary);
        cereal::BinaryInputArchive ar(s);

        ar >> out_value;

        if (rank == (int) root)
            {
            delete[] send_counts;
            delete[] displs;
            delete[] sbuf;
            }
        delete[] rbuf;
        }

    //! Wrapper around MPI_Gatherv
    template<typename T>
    void gather_v(const T& in_value, std::vector<T> & out_values, unsigned int root, const MPI_Comm mpi_comm)
        {
        int rank;
        int size;
        MPI_Comm_rank(mpi_comm, &rank);
        MPI_Comm_size(mpi_comm, &size);

        // serialize in_value
        std::stringstream s(std::ios_base::out | std::ios_base::binary);
        cereal::BinaryOutputArchive ar(s);

        ar << in_value;
        s.flush();

        // copy into send buffer
        std::string str = s.str();
        unsigned int send_count = str.length();

        int *recv_counts = NULL;
        int *displs = NULL;
        if (rank == (int) root)
            {
            out_values.resize(size);
            recv_counts = new int[size];
            displs = new int[size];
            }

        // gather lengths of buffers
        MPI_Gather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, root, mpi_comm);

        char *rbuf = NULL;
        if (rank == (int) root)
            {
            unsigned int len = 0;
            for (unsigned int i = 0; i < (unsigned int) size; i++)
                {
                displs[i] = (i > 0) ? displs[i-1] + recv_counts[i-1] : 0;
                len += recv_counts[i];
                }
            rbuf = new char[len];
            }

        // now gather actual objects
        MPI_Gatherv((void *)str.data(), send_count, MPI_BYTE, rbuf, recv_counts, displs, MPI_BYTE, root, mpi_comm);

        // on root processor, de-serialize data
        if (rank == (int) root)
            {
            for (unsigned int i = 0; i < out_values.size(); i++)
                {
                std::stringstream s(std::string(rbuf + displs[i], recv_counts[i]), std::ios_base::in | std::ios_base::binary);
                cereal::BinaryInputArchive ar(s);

                ar >> out_values[i];
                }

            delete[] displs;
            delete[] recv_counts;
            delete[] rbuf;
            }
        }

    //! Wrapper around MPI_Allgatherv
    template<typename T>
    void all_gather_v(const T& in_value, std::vector<T> & out_values, const MPI_Comm mpi_comm)
        {
        int rank;
        int size;
        MPI_Comm_rank(mpi_comm, &rank);
        MPI_Comm_size(mpi_comm, &size);

        // serialize in_value
        std::stringstream s(std::ios_base::out | std::ios_base::binary);
        cereal::BinaryOutputArchive ar(s);

        ar << in_value;
        s.flush();

        // copy into send buffer
        std::string str = s.str();
        unsigned int send_count = str.length();

        // allocate memory for buffer lengths
        out_values.resize(size);
        int *recv_counts = new int[size];
        int *displs = new int[size];

        // gather lengths of buffers
        MPI_Allgather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, mpi_comm);

        // allocate receiver buffer
        unsigned int len = 0;
        for (unsigned int i = 0; i < (unsigned int) size; i++)
            {
            displs[i] = (i > 0) ? displs[i-1] + recv_counts[i-1] : 0;
            len += recv_counts[i];
            }
        char *rbuf = new char[len];

        // now gather actual objects
        MPI_Allgatherv((void *)str.data(), send_count, MPI_BYTE, rbuf, recv_counts, displs, MPI_BYTE, mpi_comm);

        // de-serialize data
        for (unsigned int i = 0; i < out_values.size(); i++)
            {
            std::stringstream s(std::string(rbuf + displs[i], recv_counts[i]), std::ios_base::in | std::ios_base::binary);
            cereal::BinaryInputArchive ar(s);

            ar >> out_values[i];
            }

        delete[] displs;
        delete[] recv_counts;
        delete[] rbuf;
        }

    //! Wrapper around MPI_Send that handles any serializable object
    template<typename T>
    void send(const T& val,const unsigned int dest, const MPI_Comm mpi_comm)
        {
        int rank;
        MPI_Comm_rank(mpi_comm, &rank);
        if(rank == static_cast<int>(dest) ) //Quick exit, if dest is src
          return;
        char *buf = NULL;
        int recv_count;

        std::stringstream s(std::ios_base::out | std::ios_base::binary);
        cereal::BinaryOutputArchive ar(s);

        // serialize object
        ar << val;

        // do not forget to flush stream
        s.flush();

        // copy string to send buffer
        std::string str = s.str();
        recv_count = str.size();
        buf = new char[recv_count];
        str.copy(buf, recv_count);

        MPI_Send(&recv_count, 1 , MPI_INT, dest, 0, mpi_comm);

        MPI_Send(buf, recv_count, MPI_BYTE, dest, 0 , mpi_comm);

        delete[] buf;
        }

    //! Wrapper around MPI_Recv that handles any serializable object
    template<typename T>
    void recv(T& val,const unsigned int src, const MPI_Comm mpi_comm)
        {
        int rank;
        MPI_Comm_rank(mpi_comm, &rank);
        if( rank == static_cast<int>(src) ) //Quick exit if src is dest.
          return;

        int recv_count;

        MPI_Recv(&recv_count, 1, MPI_INT, src, 0 , mpi_comm, MPI_STATUS_IGNORE);

        char *buf = NULL;
        buf = new char[recv_count];

        MPI_Recv(buf, recv_count, MPI_BYTE, src, 0 , mpi_comm, MPI_STATUS_IGNORE);

        // de-serialize
        std::stringstream s(std::string(buf, recv_count), std::ios_base::in | std::ios_base::binary);
        cereal::BinaryInputArchive ar(s);
        ar >> val;

        delete[] buf;
        }
}
#endif // __MPI_INTERFACE_H__

