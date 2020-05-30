#ifndef __THERMO_PROPERTY_H__
#define __THERMO_PROPERTY_H__

#include <pyglasstools/Observable.h>
#include <Eigen/CXX11/Tensor>
#include <array>
namespace abr = Aboria;

class PYBIND11_EXPORT ThermoProperty : public Observable
{
    public:
        ThermoProperty() : Observable() {};

        ThermoProperty(std::string _name, std::string _type, bool _islocal, int _dim) 
            : Observable(_name, _type, _islocal, _dim){};
        
        virtual void accumulate(const AboriaParticles::value_type& particle_i) 
        {
        }
        
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j,
                                Eigen::Vector3d rij) 
        {
        }
        
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j,
                                Eigen::Vector3d rij, 
                                const std::shared_ptr<PairPotential>& potential)
        {
        }

        virtual void save( std::shared_ptr< MPI::LogFile > logfile, int index){}

        virtual void clear(){}
        
        virtual void divideByVolume(double vol){}
    
};

//Now, create specializations
template< int Rank, int Dim, typename AtomicObs>
class PYBIND11_EXPORT LocalProperty : public ThermoProperty
{
    private:
        AtomicObs obs;
        Eigen::Tensor<double, Rank> val;
    public:
        LocalProperty()
        {
            clear();
        };
        LocalProperty(std::string _name, std::string _type, bool _islocal) 
            : ThermoProperty(_name, _type, _islocal, Dim) 
        {
            //if (Rank == 0)
            //    val = Eigen::Tensor<double, Rank, 1>;
            //else
            //{
                //Resize so that every dimension is as long as the physical dimension
            std::array<Eigen::Index, Rank > size_array;
            size_array.fill(Dim);
            val.resize(size_array);
            //}
            clear();
        };
        
        void accumulate(const AboriaParticles::value_type& particle_i) 
        {
            obs.compute(particle_i, val);
        }
        
        void save( std::shared_ptr< MPI::LogFile > logfile, int index)
        {
            Eigen::TensorMap< Eigen::Tensor<double, 1> > tensorview(val.data(), val.size());
            logfile->write_shared(std::to_string(tensorview(index)));
        }

        void clear()
        {
            val.setZero();
        }
        
        void divideByVolume(double vol)
        {
            val = val*(1/vol);
        }
    
};

template< int Rank, int Dim, typename AtomicObs>
class PYBIND11_EXPORT ForceProperty : public ThermoProperty
{
    private:
        AtomicObs obs;
        Eigen::Tensor<double, Rank> val;
    public:
        ForceProperty()
        {
            clear();
        };
        ForceProperty(std::string _name, std::string _type, bool _islocal) 
            : ThermoProperty(_name, _type, _islocal, Dim) 
        {
            //if (Rank == 0)
            //    val = Eigen::Tensor<double, Rank, 1>;
            //else
            //{
                //Resize so that every dimension is as long as the physical dimension
            std::array<Eigen::Index, Rank > size_array;
            size_array.fill(Dim);
            val.resize(size_array);
            //}
            clear();
        };
        
        void accumulate(const AboriaParticles::value_type& particle_i, 
                        const AboriaParticles::value_type& particle_j,
                        Eigen::Vector3d rij, 
                        const std::shared_ptr<PairPotential>& potential)
        {
            obs.compute(particle_i,particle_j,rij,potential, val);
        }
        
        void save( std::shared_ptr< MPI::LogFile > logfile, int index)
        {
            /*
             * If I need this again, I'll just un-comment it
             * This is only for checking the 4-rank tensor layout!
            py::print("Full Tensor");
            for(int l = 0; l < val.dimension(3); ++l)
            {
                for(int k = 0; k < val.dimension(2); ++k)
                {
                    for(int j = 0; j < val.dimension(1); ++j)
                    {
                        for (int i = 0; i < val.dimension(0); ++i)
                        {
                            py::print(val(i,j,k,l),i,j,k,l);
                        }
                    }
                }
            }
            py::print("Flattened Tensor");
            for(int l = 0; l < val.size(); ++l)
            {
                py::print(tensorview(l),l);
            }
            */
            Eigen::TensorMap< Eigen::Tensor<double, 1> > tensorview(val.data(), val.size());
            logfile->write_shared(std::to_string(tensorview(index)));
        }

        void clear()
        {
            val.setZero();
        }
        
        void divideByVolume(double vol)
        {
            val = val*(1/vol);
        }
    
};

void export_ThermoPropertyBase(py::module& m)
{
    py::class_<ThermoProperty, Observable, std::shared_ptr< ThermoProperty > >(m,"ThermoProperty")
    ;
};

template<class T>
void export_ThermoProperty(py::module& m, const std::string& name)
{
    py::class_<T, ThermoProperty, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::string, std::string, bool>())
    .def("save",&T::save) 
    .def("clear",&T::clear) 
    ;
};

#endif //__THERMO_PROPERTY_H__
