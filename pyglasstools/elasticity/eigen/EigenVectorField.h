#ifndef __EIGEN_VECTOR_FIELD_H__
#define __EIGEN_VECTOR_FIELD_H__

#include <pyglasstools/Observable.h>
#include <pyglasstools/Manager.h>

class PYBIND11_EXPORT EigenVectorFieldBase : public Observable
{
    protected:
        std::stringstream outline;
        bool notconstructed;
    public:
        Eigen::VectorXd vectorobs;
        
        EigenVectorFieldBase() : Observable() {};
        EigenVectorFieldBase(   std::string _name, bool _islocal, int _dim)
            : Observable(_name, "VECTOR", _islocal, _dim), outline(""), notconstructed(true)
        {
        };

        ~EigenVectorFieldBase(){};
        //Virtual methods for vector fields
        virtual std::vector<double> getVectorValue(unsigned int id){return std::vector<double>();};
        virtual std::string gettostring(const int id){return outline.str();};
};

template< int Dim >
class PYBIND11_EXPORT EigenVectorField : public EigenVectorFieldBase
{
    public:
        EigenVectorField(){};
        EigenVectorField(   std::string _name, bool _islocal)
            : EigenVectorFieldBase(_name, _islocal, Dim)         
        {
        };
            //Vector Construction should be done later.
        ~EigenVectorField()
        {
        }       
        
        virtual std::vector<double> getVectorValue(unsigned int id)
        {
            //clear the stringstream
            std::vector<double> output(3);
            int idx = Dim*id;
            int idy = Dim*id+1;
            output[0] = vectorobs[idx];
            output[1] = vectorobs[idy];
            if(Dim == 3)
            {
                int idz = Dim*id+2;
                output[2] = vectorobs[idz];
            }
            return output;
        };
        
        virtual std::string gettostring(const int id)
        {
            //clear the stringstream
            outline.str( std::string() );
            outline.clear();
            int idx = Dim*id;
            int idy = Dim*id+1;
            outline << detail::to_string_sci(vectorobs[idx]) + " ";
            outline << detail::to_string_sci(vectorobs[idy]) + " ";
            
            if(Dim == 3)
            {
                int idz = Dim*id+2;
                outline << detail::to_string_sci(vectorobs[idz]) + " ";
            }
            return outline.str();
        }
};

void export_EigenVectorFieldBase(py::module& m)
{
    py::class_<EigenVectorFieldBase, std::shared_ptr<EigenVectorFieldBase> >(m,"EigenVectorFieldBase")
    .def(py::init< std::string, bool, int >())
    ;
};


template<class T>
void export_EigenVectorField(py::module& m, const std::string& name)
{
    py::class_<T, EigenVectorFieldBase, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::string, bool >())
    .def("gettostring", &T::gettostring)
    .def("getVectorValue", &T::getVectorValue)
    ;
};

#endif //__EIGEN_VECTOR_FIELD_H__
