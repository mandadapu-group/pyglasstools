#ifndef __COARSEGRAINED_FIELD_H__
#define __COARSEGRAINED_FIELD_H__

#include <pyglasstools/Observable.h>
#include <Eigen/CXX11/Tensor>
#include <array>

class PYBIND11_EXPORT CoarseGrainedField : public Observable
{
    protected:
        std::shared_ptr< std::vector< Eigen::Vector3d > > m_gridpoints;
        std::stringstream outline;
    public:
        CoarseGrainedField() : Observable() {};
        CoarseGrainedField(std::string _name, std::string _type, bool _islocal, int _dim) 
            : Observable(_name, _type, _islocal, _dim){};
        
        virtual void addGridpoints(const std::vector< Eigen::Vector3d >& gridpoints)
        {
            m_gridpoints = std::make_shared< std::vector< Eigen::Vector3d > >(gridpoints);
        };

        virtual void accumulate(const AboriaParticles::value_type& particle_i,
                                double cgval, unsigned int grid_id)
        {
        };

        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j, 
                                Eigen::Vector3d rij,
                                double bondval, unsigned int grid_id)
        {
        };

        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j, 
                                Eigen::Vector3d rij,
                                const std::shared_ptr<PairPotential>& potential,
                                double bondval, unsigned int grid_id)
        {
        };
        
        virtual void save( std::shared_ptr< MPI::LogFile > logfile, int index){}

        virtual void clear( unsigned int grid_id){}
};


//Now, create specializations
template< int Rank, int Dim, typename AtomicObs >
class PYBIND11_EXPORT LocalField : public CoarseGrainedField
{
    private:
        AtomicObs obs;
        std::vector< Eigen::Tensor<double, Rank> > val;
    public:
        LocalField(){};
        LocalField(std::string _name, std::string _type, bool _islocal, int _gridsize) 
            : CoarseGrainedField(_name, _type, _islocal, Dim)
        {
            val.resize(_gridsize);
            for (int i = 0; i < _gridsize; ++i)
                clear(i);
        };
        
        virtual void accumulate(const AboriaParticles::value_type& particle_i,
                                double cgval, unsigned int grid_id)
        {
            obs.compute_cg(particle_i, val[grid_id], cgval);
        };
        
        virtual void save( std::shared_ptr< MPI::LogFile > logfile, int index, int grid_id)
        {
            //clear the stringstream
            outline.str( std::string() );
            outline.clear();
            //for(unsigned int i = 0; i < val.size(); ++i)
            //{
            //This process might be expensive
            Eigen::TensorMap< Eigen::Tensor<double, 1> > tensorview(val[grid_id].data(), val[grid_id].size());
            
            //Input grid position 
            /*
            outline << (*m_gridpoints)[i][0] << " ";
            outline << (*m_gridpoints)[i][1] << " ";
            if (Dim == 3)
                outline << (*m_gridpoints)[i][2] << " ";
            */
            //Input the observable
            outline << tensorview(index) << " ";
            //end the input
            //outline << std::endl;
            logfile->write_shared(outline.str());
            //}
        }
        
        virtual std::string gettostring(int index, int grid_id)
        {
            //for(unsigned int i = 0; i < val.size(); ++i)
            //{
            //This process might be expensive
            //clear the stringstream
            outline.str( std::string() );
            outline.clear();
            Eigen::TensorMap< Eigen::Tensor<double, 1> > tensorview(val[grid_id].data(), val[grid_id].size());
            
            //Input grid position 
            /*
            outline << (*m_gridpoints)[i][0] << " ";
            outline << (*m_gridpoints)[i][1] << " ";
            if (Dim == 3)
                outline << (*m_gridpoints)[i][2] << " ";
            */
            //Input the observable
            outline << tensorview(index) << " ";
            //end the input
            //outline << std::endl;
            return outline.str();
            //}
        }
        
        void clear(unsigned int grid_id)
        {
            val[grid_id].setZero();
            //resize(size_array);
        } 
};

//This class has a misleading name. I should probably change it.
template< int Rank, int Dim, typename AtomicObs >
class PYBIND11_EXPORT ForceField : public CoarseGrainedField
{
    private:
        AtomicObs obs;
        std::vector< Eigen::Tensor<double, Rank> > val;
    public:
        ForceField(){};
        ForceField(std::string _name, std::string _type, bool _islocal, int _gridsize) 
            : CoarseGrainedField(_name, _type, _islocal, Dim)
        {
            val.resize(_gridsize);
            std::array<Eigen::Index, Rank > size_array;
            size_array.fill(Dim);
            for (int i = 0; i < _gridsize; ++i)
            {
                val[i].resize(size_array);
                clear(i);
            }
        };
        
        virtual void accumulate(const AboriaParticles::value_type& particle_i, 
                                const AboriaParticles::value_type& particle_j,
                                Eigen::Vector3d rij, 
                                const std::shared_ptr<PairPotential>& potential,
                                double bondval, unsigned int grid_id)
        {
            obs.compute_cg(particle_i, particle_j, rij, potential, val[grid_id], bondval);
        }
        
        virtual std::string gettostring(int index, int grid_id)
        {
            //for(unsigned int i = 0; i < val.size(); ++i)
            //{
            //This process might be expensive
            //clear the stringstream
            outline.str( std::string() );
            outline.clear();
            Eigen::TensorMap< Eigen::Tensor<double, 1> > tensorview(val[grid_id].data(), val[grid_id].size());
            
            //Input grid position 
            /*
            outline << (*m_gridpoints)[i][0] << " ";
            outline << (*m_gridpoints)[i][1] << " ";
            if (Dim == 3)
                outline << (*m_gridpoints)[i][2] << " ";
            */
            //Input the observable
            outline << tensorview(index) << " ";
            //end the input
            //outline << std::endl;
            return outline.str();
            //}
        }
        
        virtual void save( std::shared_ptr< MPI::LogFile > logfile, int index, int grid_id)
        {
            //for(unsigned int i = 0; i < val.size(); ++i)
            //{
            //This process might be expensive
            Eigen::TensorMap< Eigen::Tensor<double, 1> > tensorview(val[grid_id].data(), val[grid_id].size());
            
            //Input grid position 
            /*
            outline << (*m_gridpoints)[i][0] << " ";
            outline << (*m_gridpoints)[i][1] << " ";
            if (Dim == 3)
                outline << (*m_gridpoints)[i][2] << " ";
            */
            //Input the observable
            outline << tensorview(index) << " ";
            //end the input
            //outline << std::endl;
            logfile->write_shared(outline.str());
            
            //clear the stringstream
            outline.str( std::string() );
            outline.clear();
            //}
        }
        
        
        void clear(unsigned int grid_id)
        {
            val[grid_id].setZero();
        } 
};

void export_CoarseGrainedFieldBase(py::module& m)
{
    py::class_<CoarseGrainedField, Observable, std::shared_ptr< CoarseGrainedField > >(m,"CoarseGrainedField")
    ;
};

template<class T>
void export_CoarseGrainedField(py::module& m, const std::string& name)
{
    py::class_<T, CoarseGrainedField, std::shared_ptr<T> >(m,name.c_str())
    .def(py::init< std::string, std::string, bool, int>())
    .def("addGridpoints",&T::addGridpoints) 
    .def("gettostring", &T::gettostring)
    .def("save",&T::save) 
    .def("clear",&T::clear) 
    ;
};

#endif //__COARSEGRAINED_FIELD_H__
