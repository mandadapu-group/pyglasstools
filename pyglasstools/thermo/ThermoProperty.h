#ifndef __THERMO_PROPERTY_H__
#define __THERMO_PROPERTY_H__

#include <pyglasstools/Observable.h>
namespace abr = Aboria;

template< class Property >
class PYBIND11_EXPORT ThermoProperty : public Observable
{
    private:
        Property obs;
        Eigen::MatrixXd val;
    public:
        ThermoProperty()
        {
            clear();
        };
        ThermoProperty(std::string _name, std::string _type, bool _islocal, int _dim) 
            : Observable(_name, _type, _islocal, _dim), obs(_dim) 
        {
            clear();
        };

        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j, 
                                const std::shared_ptr<PairPotential>& potential)
        {
            val += obs.compute(particle_i, particle_j, potential);
        }

        void save( std::shared_ptr< MPI::LogFile > logfile, std::vector< unsigned int > indices)
        {
            if (type == "SCALAR")
                logfile->write_shared(std::to_string(val(0,0)));
            else if (type == "VECTOR")
                logfile->write_shared(std::to_string(val(indices[0],0)));
            else if (type == "2-TENSOR" || type == "4-TENSOR")
            {
                logfile->write_shared(std::to_string(val(indices[0],indices[1])));
            }
        } 
        void clear()
        {
            if (type == "SCALAR")
                val = Eigen::MatrixXd::Zero(1,1);
            else if (type == "VECTOR")
                val = Eigen::MatrixXd::Zero(dim,1);
            else if (type == "2-TENSOR")
                val = Eigen::MatrixXd::Zero(dim,dim);
            else if (type == "4-TENSOR")
                val = Eigen::MatrixXd::Zero((int)3*(dim-1),(int)3*(dim-1));
        }
        
        void divideByVolume(double vol)
        {
            val /= vol;
        }
    
};

template<class T>
void export_ThermoProperty(py::module& m, const std::string& name)
{
    py::class_<T, Observable, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::string, std::string, bool, int >())
    .def("save",&T::save) 
    .def("clear",&T::clear) 
    ;
};
#endif //__THERMO_PROPERTY_H__
