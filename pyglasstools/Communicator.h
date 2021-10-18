#ifndef __COMMUNICATOR_H__
#define __COMMUNICATOR_H__

#include <Eigen/Core>

#include <pybind11/pybind11.h>
#include <stdlib.h>

/*Class: Communicator
 *A Base class for communicator, so that purely serial and MPI version of pyglasstools can work in the same way
 *Default methods implemented are rank and size information 
 */

class PYBIND11_EXPORT Communicator
{
    public:
        //! Constructor with externally provided MPI communicator (only in MPI enabled builds)
        Communicator(){};
        
        //! Destructor
        virtual ~Communicator() {};
        
        //! Return the rank of this processor in the partition
        unsigned int getRank() const
        {
            return 0;
        };
        
        //! Return the number of ranks in this partition
        unsigned int getNRanks() const
        {
            return 1;
        };
        
        //! Returns true if this is the root processor
        bool isRoot() const
        {
            return true;
        };

        void abort(int error_code) const
        {
            // delay for a moment to give time for error messages to print
            int msec = 1000;
            usleep(msec*1000);
            std::exit(error_code);
        };
};

void export_Communicator(pybind11::module& m)
{
    pybind11::class_<Communicator, std::shared_ptr<Communicator> > communicator(m,"Communicator");
    communicator.def(pybind11::init< >())
    .def("getRank", &Communicator::getRank)
    .def("getNRanks", &Communicator::getNRanks)
    .def("isRoot", &Communicator::isRoot)
    .def("abort", &Communicator::abort)
    ;
}

#endif
