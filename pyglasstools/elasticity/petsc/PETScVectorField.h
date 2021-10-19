#ifndef __PETSC_VECTOR_FIELD_H__
#define __PETSC_VECTOR_FIELD_H__

#include <pyglasstools/Observable.h>
#include <pyglasstools/Manager.h>
#include <pyglasstools/MPICommunicator.h>
#include <petscmat.h>

class PYBIND11_EXPORT PETScVectorFieldBase : public Observable
{
    protected:
        std::shared_ptr< MPI::ParallelCommunicator > m_comm;
        std::stringstream outline;
        bool notconstructed;
    public:
        Vec vectorobs;
        
        PETScVectorFieldBase() : Observable() {};
        PETScVectorFieldBase(   std::string _name, bool _islocal, int _dim, 
                                std::shared_ptr< MPI::ParallelCommunicator > comm) 
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
        //PetscInt Istart;
        //PetscInt Iend;
        std::vector<double> iovector;

        //std::vector< PetscInt > globalIstart; 
        //std::vector< PetscInt > globalIend;
        
        PETScVectorField(){};
        PETScVectorField(   std::string _name, bool _islocal, 
                            std::shared_ptr< MPI::ParallelCommunicator > comm) 
            : PETScVectorFieldBase(_name, _islocal, Dim,comm)         
        {
        };
            //Vector Construction should be done later.
        ~PETScVectorField()
        {
            if (!notconstructed)
                VecDestroy(&vectorobs);
        }       
        void createVector(Mat A)
        {
            if (notconstructed)
            {
                MatCreateVecs(A,NULL,&vectorobs);
                notconstructed =  false;
            }
        }
        
        void scatterVector(int firstparticleid, int localsize)
        {
            VecScatter  scatter;      /* scatter context */
            PetscScalar *values;       //a temporary array for which I can store inside iovector
            PetscInt    *idx_from = new PetscInt[2*localsize];
            PetscInt    *idx_to = new PetscInt[2*localsize];
            for (int i = 0; i < 2*localsize; ++i)
            {
                idx_from[i] = 2*firstparticleid+i;
                idx_to[i] = i;
                //py::print(idx_from[i],"for Process ",m_comm->getRank(), "with first particle id ",firstparticleid, "and local size ",2*localsize);
            } 
            //Create temporary sequential vector
            Vec x;
            VecCreateSeq(PETSC_COMM_SELF,2*localsize,&x);
            
            IS          from, to;     /* index sets that define the scatter */
            ISCreateGeneral(PETSC_COMM_SELF,2*localsize,idx_from,PETSC_COPY_VALUES,&from);
            ISCreateGeneral(PETSC_COMM_SELF,2*localsize,idx_to,PETSC_COPY_VALUES,&to);
            
            /* Now scatter our array! */
            VecScatterCreate(vectorobs,from,x,to,&scatter);
            VecScatterBegin(scatter,vectorobs,x,INSERT_VALUES,SCATTER_FORWARD);
            VecScatterEnd(scatter,vectorobs,x,INSERT_VALUES,SCATTER_FORWARD);
            VecGetArray(x,&values);
            iovector.resize(2*localsize); 
            for (int i = 0; i < 2*localsize; ++i)
            {
                iovector[i] = values[i];
            }
            /* Finally, destroy these values */
            ISDestroy(&from);
            ISDestroy(&to);
            VecScatterDestroy(&scatter);
            VecDestroy(&x);
            delete [] idx_from;
            delete [] idx_to;
        } 
        
        Vec getVector()
        {
            return vectorobs;
        }
        /*
        int whichRankOwnsThis(int& id)
        {
            for(int i = 0; i < m_comm->getSizeGlobal(); i++)
            {
                if (id >= globalIstart[i] && id < globalIend[i])
                {
                    return i;
                }
            }
            return 0;
        }
        */
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
                //double valx, valy;
                outline << detail::to_string_sci(iovector[idx]) + " ";
                outline << detail::to_string_sci(iovector[idy]) + " ";
                //}
                //VecGetValues(vectorobs,1,&idy,&valy);
                //outline << detail::to_string_sci(valy) + " ";
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
    .def(py::init< std::string, bool, int, std::shared_ptr< MPI::ParallelCommunicator > >())
    ;
};


template<class T>
void export_PETScVectorField(py::module& m, const std::string& name)
{
    py::class_<T, PETScVectorFieldBase, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::string, bool, std::shared_ptr< MPI::ParallelCommunicator > >())
    .def("gettostring", &T::gettostring)
    .def("getVectorValue", &T::getVectorValue)
    .def("scatterVector", &T::scatterVector)
    ;
};

#endif //__PETSC_VECTOR_FIELD_H__
