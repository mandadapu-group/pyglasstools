// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.
#ifndef __MPI_INTERFACE_H__
#define __MPI_INTERFACE_H__

// ensure that HOOMDMath.h is the first thing included
#include <mpi.h>
#include <sstream>
#include <vector>

/*! \file Communicator.h
    \brief Declares Communicator, which initializes the MPI environment
*/
#include <Eigen/Core>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp> 
#include <cereal/archives/binary.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

//for serialization of eigen types
namespace cereal
{
  template <class Archive, class Derived> inline
    typename std::enable_if<traits::is_output_serializable<BinaryData<typename Derived::Scalar>, Archive>::value, void>::type
    save(Archive & ar, Eigen::PlainObjectBase<Derived> const & m){
      typedef Eigen::PlainObjectBase<Derived> ArrT;
      if(ArrT::RowsAtCompileTime==Eigen::Dynamic) ar(m.rows());
      if(ArrT::ColsAtCompileTime==Eigen::Dynamic) ar(m.cols());
      ar(binary_data(m.data(),m.size()*sizeof(typename Derived::Scalar)));
    }

  template <class Archive, class Derived> inline
    typename std::enable_if<traits::is_input_serializable<BinaryData<typename Derived::Scalar>, Archive>::value, void>::type
    load(Archive & ar, Eigen::PlainObjectBase<Derived> & m){
      typedef Eigen::PlainObjectBase<Derived> ArrT;
      Eigen::Index rows=ArrT::RowsAtCompileTime, cols=ArrT::ColsAtCompileTime;
      if(rows==Eigen::Dynamic) ar(rows);
      if(cols==Eigen::Dynamic) ar(cols);
      m.resize(rows,cols);
      ar(binary_data(m.data(),static_cast<std::size_t>(rows*cols*sizeof(typename Derived::Scalar))));
    }
}

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
            MPI_Comm getWorldCommunicator() const
            {
                return m_comm_world;
            }

            void mpiabort(int error_code)
            {
                if( getSizeGlobal() > 1)
                {
                    // delay for a moment to give time for error messages to print
                    int msec = 1000;
                    usleep(msec*1000);
                    MPI_Abort(m_comm, error_code);
                }
            }

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
            unsigned int getSizeGlobal() const
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
                return getSizeGlobal()/m_n_rank;
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
            
            //! Abort MPI runs
            void abort(unsigned int milliseconds)
            {
                if(getSizeGlobal() > 1)
                {
                    // delay for a moment to give time for error messages to print
                    usleep(1000 * milliseconds);
                    MPI_Abort(m_comm_world, MPI_ERR_OTHER);
                }
            }
            //! Wrapper around MPI_Bcast that handles any serializable object
            template<typename T>
            void bcast(T& val, unsigned int root)
                {
                int rank;
                MPI_Comm_rank(m_comm, &rank);

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

                MPI_Bcast(&recv_count, 1, MPI_INT, root, m_comm);
                if (rank != (int) root)
                    buf = new char[recv_count];

                MPI_Bcast(buf, recv_count, MPI_BYTE, root, m_comm);

                if (rank != (int)root)
                    {
                    // de-serialize
                    std::stringstream s(std::string(buf, recv_count), std::ios_base::in | std::ios_base::binary);
                    cereal::BinaryInputArchive ar(s);

                    ar >> val;
                    }

                delete[] buf;
                }

            template<typename T>
            T pybcast(T& val, unsigned int root)
            {
                T new_val = val;
                bcast(new_val,root);
                return new_val;
            }
            //! Wrapper around MPI_Scatterv that scatters a vector of serializable objects
            template<typename T>
            void scatter_v(const std::vector<T>& in_values, T& out_value, unsigned int root)
                {
                int rank;
                int size;
                MPI_Comm_rank(m_comm, &rank);
                MPI_Comm_size(m_comm, &size);

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
                MPI_Scatter(send_counts, 1, MPI_INT, &recv_count, 1, MPI_INT, root, m_comm);

                // allocate receive buffer
                char *rbuf = new char[recv_count];

                // scatter actual data
                MPI_Scatterv(sbuf, send_counts, displs, MPI_BYTE, rbuf, recv_count, MPI_BYTE, root, m_comm);

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

            template<typename T>
            T pyscatter_v(const std::vector<T>& in_values, unsigned int root)
            {
                T out_value;
                scatter_v(in_values,out_value,root);
                return out_value;
            }
            //! Wrapper around MPI_Gatherv
            template<typename T>
            void gather_v(const T& in_value, std::vector<T> & out_values, unsigned int root)
                {
                int rank;
                int size;
                MPI_Comm_rank(m_comm, &rank);
                MPI_Comm_size(m_comm, &size);

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
                MPI_Gather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, root, m_comm);

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
                MPI_Gatherv((void *)str.data(), send_count, MPI_BYTE, rbuf, recv_counts, displs, MPI_BYTE, root, m_comm);

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
            void all_gather_v(const T& in_value, std::vector<T> & out_values)
                {
                int rank;
                int size;
                MPI_Comm_rank(m_comm, &rank);
                MPI_Comm_size(m_comm, &size);

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
                MPI_Allgather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, m_comm);

                // allocate receiver buffer
                unsigned int len = 0;
                for (unsigned int i = 0; i < (unsigned int) size; i++)
                    {
                    displs[i] = (i > 0) ? displs[i-1] + recv_counts[i-1] : 0;
                    len += recv_counts[i];
                    }
                char *rbuf = new char[len];

                // now gather actual objects
                MPI_Allgatherv((void *)str.data(), send_count, MPI_BYTE, rbuf, recv_counts, displs, MPI_BYTE, m_comm);

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
            
            //! Wrapper around MPI_Allgatherv
            template<typename T>
            std::vector<T> pyall_gather_v(const T& in_value)
            {
                std::vector<T> out_values;
                all_gather_v(in_value, out_values);
                return out_values;
            }
            //! Wrapper around MPI_Send that handles any serializable object
            template<typename T>
            void send(const T& val,const unsigned int dest, const unsigned int tag)
                {
                int rank;
                MPI_Comm_rank(m_comm, &rank);
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

                MPI_Send(&recv_count, 1 , MPI_INT, dest, tag, m_comm);

                MPI_Send(buf, recv_count, MPI_BYTE, dest, tag , m_comm);

                delete[] buf;
                }

            //! Wrapper around MPI_Recv that handles any serializable object
            template<typename T>
            void recv(T& val,const unsigned int src, const unsigned int tag)
                {
                int rank;
                MPI_Comm_rank(m_comm, &rank);
                if( rank == static_cast<int>(src) ) //Quick exit if src is dest.
                  return;

                int recv_count;

                MPI_Recv(&recv_count, 1, MPI_INT, src, tag , m_comm, MPI_STATUS_IGNORE);

                char *buf = NULL;
                buf = new char[recv_count];

                MPI_Recv(buf, recv_count, MPI_BYTE, src, tag , m_comm, MPI_STATUS_IGNORE);

                // de-serialize
                std::stringstream s(std::string(buf, recv_count), std::ios_base::in | std::ios_base::binary);
                cereal::BinaryInputArchive ar(s);
                ar >> val;

                delete[] buf;
                }
            
            template<typename T>
            T pyrecv(const unsigned int src, const unsigned int tag)
                {
                    T val;
                    recv(val,src,tag);
                    return val;
                }

        protected:
            MPI_Comm m_comm;                    //!< The MPI communicator
            MPI_Comm m_comm_world;              //!< The world communicator
            
            unsigned int m_rank;                //!< Rank of this processor (0 if running in single-processor mode)
            unsigned int m_n_rank;              //!< Ranks per partition
        };
}

//! Exports Communicator to python
void export_MPICommunicator(pybind11::module& m);

#endif
