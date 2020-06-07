#ifndef __PETSC_GLOBAL_PROPERTY_H__
#define __PETSC_GLOBAL_PROPERTY_H__

#include <pyglasstools/Observable.h>
#include <pyglasstools/Manager.h>
#include <array>
#include <Eigen/CXX11/Tensor>

class PYBIND11_EXPORT PETScGlobalPropertyBase : public Observable
{
    protected:
        std::shared_ptr< MPI::Communicator > m_comm;
        std::stringstream outline;
    public:
        PETScGlobalPropertyBase() : Observable() {};
        PETScGlobalPropertyBase(    std::string _name, std::string _type, bool _islocal, int _dim, 
                            std::shared_ptr< MPI::Communicator > comm) 
            : Observable(_name, _type, _islocal, _dim), m_comm(comm), outline("")
        {
        };

        ~PETScGlobalPropertyBase(){};

        virtual void setValue(double valid, const int id = 0){}
        virtual void addValue(double valid, int id){}
        virtual std::vector<double> getValue(){return std::vector<double>(1);}
        //Virtual methods for vector fields
        virtual std::string gettostring(int id = 0){outline.str(""); return outline.str();};
        virtual void save( std::shared_ptr< MPI::LogFile > logfile, int index = 0){}
        virtual int getDimension(int dim){ return 0;};
        virtual void multiplyBy(double scalar){}
};


class PYBIND11_EXPORT PETScGlobalScalar : public PETScGlobalPropertyBase
{
    protected:
        double val;
    public:
        
        PETScGlobalScalar(){};
        PETScGlobalScalar(    std::string _name, std::string _type, bool _islocal, 
                                std::shared_ptr< MPI::Communicator > comm) 
            : PETScGlobalPropertyBase(_name, _type, _islocal, 2 ,comm), val(0)
        {
        };

        ~PETScGlobalScalar(){};
        
        //Virtual methods for vector fields
        virtual void setValue(double valid, const int id = 0)
        {
            val = valid;
        }
        virtual std::string gettostring()
        {
            return std::to_string(val);
        } 
        void save( std::shared_ptr< MPI::LogFile > logfile)
        {
            logfile->write_shared(std::to_string(val));
        }
};

template< int Rank, int Dim >
class PYBIND11_EXPORT PETScGlobalProperty : public PETScGlobalPropertyBase
{
    protected:
        Eigen::Tensor<double, Rank> val;
    public:
        
        PETScGlobalProperty(){};
        PETScGlobalProperty(    std::string _name, std::string _type, bool _islocal, 
                                std::shared_ptr< MPI::Communicator > comm) 
            : PETScGlobalPropertyBase(_name, _type, _islocal, Dim,comm)
        {
            std::array<Eigen::Index, Rank > size_array;
            size_array.fill(Dim);
            val.resize(size_array);
            clear();
        };

        ~PETScGlobalProperty(){};

        virtual void setValue(double valid, int id)
        {
            double* val_data = val.data();
            val_data[id] = valid;
        }
        
        virtual void addValue(double valid, int id)
        {
            double* val_data = val.data();
            val_data[id] += valid;
        }
        
        virtual void multiplyBy(double scalar)
        {
            val = val*scalar;
        }
        
        virtual std::vector<double> getValue()
        {
            std::vector<double> outval(val.data(),val.data()+val.size());
            return outval;
        }
        //Virtual methods for vector fields
        virtual std::string gettostring(const int id)
        {
            Eigen::TensorMap< Eigen::Tensor<double, 1> > tensorview(val.data(), val.size());
            return std::to_string(tensorview(id));
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
        virtual int getDimension(int dim)
        { 
            return val.dimension(dim);
        };
};


void export_PETScGlobalPropertyBase(py::module& m)
{
    py::class_<PETScGlobalPropertyBase, std::shared_ptr<PETScGlobalPropertyBase> >(m,"PETScGlobalPropertyBase")
    .def(py::init< std::string, std::string, bool, int, std::shared_ptr< MPI::Communicator > >())
    ;
};

template<class T>
void export_PETScGlobalProperty(py::module& m, const std::string& name)
{
    py::class_<T, PETScGlobalPropertyBase, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::string, std::string, bool, std::shared_ptr< MPI::Communicator > >())
    .def("setValue", &T::setValue)
    .def("getValue", &T::getValue)
    .def("gettostring",&T::gettostring)    
    .def("save",&T::save)
    ;
};

#endif //__PETSC_GLOBAL_PROPERTY_H__
