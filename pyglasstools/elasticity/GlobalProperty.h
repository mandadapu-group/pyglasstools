#ifndef __GLOBAL_PROPERTY_H__
#define __GLOBAL_PROPERTY_H__

#include <pyglasstools/Observable.h>
#include <pyglasstools/Manager.h>
#include <array>
#include <Eigen/CXX11/Tensor>

/*
 * Base Class for storing and accessing global observables, i.e., not vector fields.
 * 
 * The class gives a set of virtual methods for accessing and saving observables
 * that must be followed for any derived classes.
 */
class PYBIND11_EXPORT GlobalPropertyBase : public Observable
{
    protected:
        std::stringstream outline;
    public:
        GlobalPropertyBase() : Observable() {};
        
        /*
         * Parametrized constructor. List of arguments
         * _name: human-readable name of the observable, e.g., "Pressure" or "Stress"
         * _type: human-readable name of the mathematical type of the observable, e.g., "Vector" or "2nd Rank Tensor"
         * _islocal: boolean addressing whether an observable requires neighboring particles for computations or not. _islocal=false implies pair computations are needed.  
         * _dim: physical dimension of the observable. Supports only 2 or 3.
         * comm: a shared pointer to the MPI communicator.
        */
        GlobalPropertyBase( std::string _name, std::string _type, bool _islocal, int _dim)
            : Observable(_name, _type, _islocal, _dim), outline("")
        {
        };

        ~GlobalPropertyBase(){};
        
        /*
         * Sets the value of the observable. Args:
         * valid: the value of a component of the observable
         * id: index of the component. 
        */
        virtual void setValue(double valid, const int id = 0){};
        
        /*
         * Adds to an existing component an additional value. Args:
         * valid: the value to be added to an existing component of the observable
         * id: index of the component. 
        */
        virtual void addValue(double valid, int id){};

        /*
         * Get the entire observable as a double vecttor.
         * All observables, whether it's a scalar, vector, or a tensor, must be returned by this type!
        */
        virtual std::vector<double> getValue()
        {
            return std::vector<double>(1);
        };
        
        /*
         * Gets a component of an observable and return as a string.
         * This is mostly useful when printing and doing print debugging. Args:
         * id: index of the component
        */
        virtual std::string gettostring(int id = 0)
        {
            outline.str(""); 
            return outline.str();
        };
        
        /* 
         * Save a component on the observable to a logfile. Args: 
         * logfile: a class managing I/O of a file. See pyglasstools/MPIFile.h
         * index: index of the component.
         */
        virtual void save( std::shared_ptr< BaseLogFile > logfile, int index = 0){};
        
        /* 
         * Return the dimensionality and size of the observable, i.e., whether it's a higher rank tensor, etc. Args:
         * dim: physical dimension we are operating on. 
         */
        virtual int getDimension(int dim)
        {
            return 0;
        };
        
        /*
         * Multiply the entire observable by some scalar value. Args: 
         * scalar: the scalar to multiply to
         */
        virtual void multiplyBy(double scalar){};
};


/*
 * Derived Class for storing and accessing scalars
 *
 * Note: In all our codes so far, it suffices to now how to set values, get them as strings,
 * and save it in a logile. 
 */
class PYBIND11_EXPORT GlobalScalar : public GlobalPropertyBase
{
    protected:
        //The stored scalar. 
        double val;

    public:
        
        GlobalScalar(){};
        GlobalScalar(std::string _name, std::string _type, bool _islocal, int _dim)//, std::shared_ptr< MPI::Communicator > comm) 
            : GlobalPropertyBase(_name, _type, _islocal, _dim), val(0)
        {
        };

        ~GlobalScalar(){};
        
        void setValue(double valid, const int id = 0)
        {
            val = valid;
        }
        std::string gettostring()
        {
            return std::to_string(val);
        } 
        void save( std::shared_ptr< BaseLogFile > logfile)
        {
            logfile->write(detail::to_string_sci(val));
        }
};

/*
class PYBIND11_EXPORT VectorField : public GlobalPropertyBase
{
    protected:
        Eigen::VectorXd val;
    public:
        
        VectorField(){};
        VectorField(    std::string _name, std::string _type, bool _islocal, int _dim, 
                                std::shared_ptr< MPI::Communicator > comm) 
            : GlobalPropertyBase(_name, _type, _islocal, _dim, comm)
        {
            val.resize(2);
            clear();
        };

        ~VectorField(){};

        virtual void setValue(double valid, int id)
        {
            val[id] = valid;
        }
        
        virtual void addValue(double valid, int id)
        {
            val[id] += valid;
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
            return std::to_string(val[id]);
        } 
        void save( std::shared_ptr< BaseLogFile > logfile, int index)
        {
            logfile->write(std::to_string(val[index]));
        }
        void clear()
        {
            val.setZero();
        }
};
*/


/*
 * Derived template class for storing and accessing vectorial and higher-rank tensorial observables.
 *
 * Note: we use the Eigen::Tensor class for Access and read operations are typically done by flattening the tensor first.
 */
template< int Rank, int Dim >
class PYBIND11_EXPORT GlobalProperty : public GlobalPropertyBase
{
    protected:

        //The stored tensor
        Eigen::Tensor<double, Rank> val;
    public:
        
        GlobalProperty(){};
        GlobalProperty(std::string _name, std::string _type, bool _islocal)//, std::shared_ptr< MPI::Communicator > comm) 
            : GlobalPropertyBase(_name, _type, _islocal, Dim)//, comm)
        {
            std::array<Eigen::Index, Rank > size_array;
            size_array.fill(Dim);
            val.resize(size_array);
            clear();
        };

        ~GlobalProperty(){};

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
        
        virtual std::string gettostring(const int id)
        {
            Eigen::TensorMap< Eigen::Tensor<double, 1> > tensorview(val.data(), val.size());
            return std::to_string(tensorview(id));
        }; 
        
        void save( std::shared_ptr< BaseLogFile > logfile, int index)
        {
            Eigen::TensorMap< Eigen::Tensor<double, 1> > tensorview(val.data(), val.size());
            logfile->write(std::to_string(tensorview(index)));
        };
        

        virtual int getDimension(int dim)
        { 
            return val.dimension(dim);
        };

        /* 
         * Clear all observables by setting all values to zero
         */
        void clear()
        {
            val.setZero();
        };
};

/*
 * Helper function to export GlobalPropertyBase class to Python
 */
void export_GlobalPropertyBase(py::module& m)
{
    py::class_<GlobalPropertyBase, std::shared_ptr<GlobalPropertyBase> >(m,"GlobalPropertyBase")
    .def(py::init< std::string, std::string, bool, int >())
    ;
};

/*
 * Helper function to export GlobalScalar class to Python.
 * Only select methods are exposed to the Python side. 
 */
void export_GlobalScalar(py::module& m, const std::string& name)
{
    py::class_<GlobalScalar, GlobalPropertyBase, std::shared_ptr<GlobalScalar> >(m,name.c_str())
    .def(py::init< std::string, std::string, bool, int >())
    .def("setValue", &GlobalScalar::setValue)
    .def("getValue", &GlobalScalar::getValue)
    .def("gettostring",&GlobalScalar::gettostring)    
    .def("save",&GlobalScalar::save)
    ;
};

/*
 * Helper function to export derived GlobalProperty classes to Python.
 * The T argument refers to specialized template class that we will export.
 * Only select methods are exposed to the Python side. 
 */
template<class T>
void export_GlobalProperty(py::module& m, const std::string& name)
{
    py::class_<T, GlobalPropertyBase, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::string, std::string, bool >())
    .def("setValue", &T::setValue)
    .def("getValue", &T::getValue)
    .def("gettostring",&T::gettostring)    
    .def("save",&T::save)
    ;
};

#endif //__GLOBAL_PROPERTY_H__
