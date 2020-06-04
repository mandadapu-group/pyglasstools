#ifndef __COARSEGRAINED_FIELD_H__
#define __COARSEGRAINED_FIELD_H__

#include <pyglasstools/Observable.h>
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

        virtual void setValue(double valid, const int id){}
        virtual std::vector<double> getValue(){return std::vector<double>(1);}
        //Virtual methods for vector fields
        virtual std::string gettostring(const int id){return outline.str();};
        virtual void save( std::shared_ptr< MPI::LogFile > logfile, int index){}
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
};

class PYBIND11_EXPORT PETScVectorFieldBase : public Observable
{
    protected:
        std::shared_ptr< MPI::Communicator > m_comm;
        std::stringstream outline;
        bool notconstructed;
    public:
        Vec vectorobs;
        
        PETScVectorFieldBase() : Observable() {};
        PETScVectorFieldBase(   std::string _name, bool _islocal, int _dim, 
                                std::shared_ptr< MPI::Communicator > comm) 
            : Observable(_name, "VECTOR", _islocal, _dim), m_comm(comm), outline(""), notconstructed(true)
        {
        };

        ~PETScVectorFieldBase(){};
        //Virtual methods for vector fields
        virtual void createVector(Mat A){} 
        virtual void copyVector(Vec x){}; 
        virtual std::vector<double> getVectorValue(unsigned int id){return std::vector<double>();};
        virtual std::string gettostring(const int id){return outline.str();};
};

template< int Dim >
class PYBIND11_EXPORT PETScVectorField : public PETScVectorFieldBase
{
    public:
        PETScVectorField(){};
        PETScVectorField(   std::string _name, bool _islocal, 
                            std::shared_ptr< MPI::Communicator > comm) 
            : PETScVectorFieldBase(_name, _islocal, Dim,comm)         
        {
        };
            //Vector Construction should be done later.
        ~PETScVectorField()
        {
            if (!notconstructed)
                VecDestroy(&vectorobs);
        }       
        //void copyVector(Mat A, Vec x)
        
        void createVector(Mat A)
        {
            if (notconstructed)
            {
                MatCreateVecs(A,NULL,&vectorobs);
                notconstructed =  false;
            }
        }
        Vec getVector()
        {
            return vectorobs;
        }
        virtual std::vector<double> getVectorValue(unsigned int id)
        {
            //clear the stringstream
            std::vector<double> output(3);
            if (!notconstructed)
            {
                
                PetscInt idx = Dim*id;
                PetscInt idy = Dim*id+1;
                double valx, valy;
                VecGetValues(vectorobs,1,&idx,&valx);
                output[0] = valx;
                VecGetValues(vectorobs,1,&idy,&valy);
                output[1] = valy;
                if(Dim == 3)
                {
                    double valz;
                    PetscInt idz = Dim*id+2;
                    VecGetValues(vectorobs,1,&idz,&valz);
                    output[2] = valz;
                }
                return output;
            }
            else
            {
                //SHoudl do a warning
                return output;
            }
        };
        
        virtual std::string gettostring(const int id)
        {
            //clear the stringstream
            outline.str( std::string() );
            outline.clear();
            if (!notconstructed)
            {
                
                PetscInt idx = Dim*id;
                PetscInt idy = Dim*id+1;
                double valx, valy;
                VecGetValues(vectorobs,1,&idx,&valx);
                outline << detail::to_string_sci(valx) + " ";
                VecGetValues(vectorobs,1,&idy,&valy);
                outline << detail::to_string_sci(valy) + " ";
                if(Dim == 3)
                {
                    double valz;
                    PetscInt idz = Dim*id+2;
                    VecGetValues(vectorobs,1,&idz,&valz);
                    outline << detail::to_string_sci(valz) + " ";
                }
                return outline.str();
            }
            else
            {
                return outline.str();
            }
        }
};

void export_PETScVectorFieldBase(py::module& m)
{
    py::class_<PETScVectorFieldBase, std::shared_ptr<PETScVectorFieldBase> >(m,"PETScVectorFieldBase")
    .def(py::init< std::string, bool, int, std::shared_ptr< MPI::Communicator > >())
    ;
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

template<class T>
void export_PETScVectorField(py::module& m, const std::string& name)
{
    py::class_<T, PETScVectorFieldBase, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::string, bool, std::shared_ptr< MPI::Communicator > >())
    .def("gettostring", &T::gettostring)
    .def("getVectorValue", &T::getVectorValue)
    ;
};

#endif //__COARSEGRAINED_FIELD_H__
