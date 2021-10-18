#ifndef __MPI_COMMUNICATOR_H__
#define __MPI_COMMUNICATOR_H__

// ensure that HOOMDMath.h is the first thing included
#include <mpl/mpl.hpp>
#include <sstream>
#include <vector>

/*
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp> 
#include <cereal/archives/binary.hpp>
*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>



namespace MPI
{
    /*Class: ParallelCommunicator
     * A light wrapper to mpl's Communnicator class
     * 
     */
    class PYBIND11_EXPORT ParallelCommunicator : public Communicator
    {
        protected:
            mpl::communicator m_comm; 
        public:
            //! Constructor with externally provided MPI communicator (only in MPI enabled builds)
            /*! \param mpi_comm world MPI communicator
             */
            ParallelCommunicator() : m_comm(mpl::environment::comm_world());
            
            //! Destructor
            virtual ~ParallelCommunicator() {};

            //! Returns the MPI communicator
            MPI_Comm getCommunicator() const;

            //! Return the rank of this processor in the partition
            unsigned int getRank() const;
            
            //! Return the number of ranks in this partition
            unsigned int getNRanks() const;
            
            //! Returns true if this is the root processor
            bool isRoot() const;

            //! Perform a job-wide MPI barrier
            void barrier();

            //! Perform a job-wide MPI Abort
            void abort(int error_code);
            
            //! Wrapper around MPI_Bcast that handles any serializable object
            template<typename T>
            inline void bcast(T& val, unsigned int root)
            {
                m_comm.bcast(root,val.data());
            }

            //! Wrapper around bcast that works from the Python side
            template<typename T>
            inline T pybcast(T& val, unsigned int root)
            {
                T new_val = val;
                this->bcast(new_val,root);
                return new_val;
            }
            
            //! Wrapper around MPI_Scatterv that scatters a vector of serializable objects
            template<typename T>
            inline void scatter_v(const std::vector<T>& in_values, T& out_value, unsigned int root)
            {
                m_comm.scatter(root, in_values.data(), out_value); 
            }

            //! Wrapper around scatter_v that works from the Python side
            template<typename T>
            T pyscatter_v(const std::vector<T>& in_values, unsigned int root)
            {
                T out_value;
                this->scatter_v(in_values,out_value,root);
                return out_value;
            }
            
            //! Wrapper around MPI_Gatherv
            template<typename T>
            inline void gather_v(const T& in_value, std::vector<T> & out_values, unsigned int root)
            {
                m_comm.gather(root,in_value, out_values.data());
            }

            //! Wrapper around MPI_Allgatherv
            template<typename T>
            inline void all_gather_v(const T& in_value, std::vector<T> & out_values)
            {
                m_comm.allgather(in_value, out_values.data());
            }
            
            //! Wrapper around all_gatherv that works from the Python side
            template<typename T>
            inline std::vector<T> pyall_gather_v(const T& in_value)
            {
                std::vector<T> out_values;
                this->all_gather_v(in_value, out_values);
                return out_values;
            }

            //! Wrapper around MPI_Send that handles any serializable object
            template<typename T>
            inline void send(const T& val,const unsigned int dest, const unsigned int tag)
            {
                m_comm.send(val,dest,tag);
            }

            //! Wrapper around MPI_Recv that handles any serializable object
            template<typename T>
            inline void recv(T& val,const unsigned int src, const unsigned int tag)
            {
                m_comm.recv(val,src,tag);
            }
            
            //!Wrapper around recv that works from the Python side
            template<typename T>
            inline T pyrecv(const unsigned int src, const unsigned int tag)
            {
                T val;
                this->recv(val,src,tag);
                return val;
            }

    };
}

//! Exports Communicator to python
void export_MPICommunicator(pybind11::module& m);

#endif
