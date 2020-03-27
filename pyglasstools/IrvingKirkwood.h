#include <vector>
#include <cmath>

//#include "PairPotential.h"
#include <Aboria.h>
#include <algorithm>
#include <memory>
#include "MathSTLVectors.h"
#include "ParticleSystem.h"
#include "CoarseGrainFunction.h"
#include "SimBox.h"

#include "extern/pybind11/include/pybind11/pybind11.h"
#include "extern/pybind11/include/pybind11/stl.h"
#include "extern/pybind11/include/pybind11/eigen.h"
namespace py = pybind11;
#include <Eigen/Dense>
#include <Eigen/StdVector>
//#include <omp.h>
using namespace Eigen;

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector3d);

class PYBIND11_EXPORT GridPoints
{
    public:
        GridPoints() : totsize(0){};
        GridPoints(const std::vector< Vector3d >& _points)
            : totsize(_points.size()), points(_points)
        {};
        ~GridPoints(){};
        
        unsigned int totsize;
        std::vector< Vector3d > points;
};

void export_GridPoints(py::module& m)
{
    py::class_<GridPoints, std::shared_ptr<GridPoints> >(m,"GridPoints")
    .def(py::init<std::vector< Vector3d > >())
    ;
};
        
class PYBIND11_EXPORT IrvingKirkwood
{
    public:
        IrvingKirkwood( const std::vector< Vector3d >& gridpoints, 
                        std::shared_ptr< ParticleSystem > sysdata, 
                        std::shared_ptr< CoarseGrainFunction > cgfunc)
                        :   grid(gridpoints), m_sysdata(sysdata), m_cgfunc(cgfunc)
        {
            //Resize the coarse grained fields according to the total number of grid points
            rho.resize(grid.totsize);
            Tv.resize(grid.totsize);
            Tk.resize(grid.totsize);
           
            //py::print("What is the value of rcut", m_sysdata->potential->getScaledRcut()); 
            //Define the largest possible radius cut-off 
            force_rcut = m_sysdata->potential->getScaledRcut();
            //py::print("The valur for scaled rcut is", force_rcut);
            maxforce_rcut = force_rcut;
            maxforce_rcut *= *std::max_element(  std::begin(get<diameter>(m_sysdata->particles)), 
                                                std::end(get<diameter>(m_sysdata->particles)) );
            //Define the fixed-length cut-off radius for the coarse graining function
            //py::print("Maximum for rcut is", maxforce_rcut);
            cg_rcut = m_cgfunc->getRcut();
            //py::print("Stored coarse graining function cut", cg_rcut);

            //Check if velocities are supplied or not
            hasvelocities = sysdata->haveVelocities();
            totTvxy = 0;
        };
        ~IrvingKirkwood(){};
        virtual void compute();
        virtual void globalcompute();
        std::vector<double> getRho()
        {
            return rho;
        };
        std::vector<double> getTv_xy()
        {
            return Tv.XY;
        };
        double getTvxy()
        {
            return totTvxy/m_sysdata->simbox->getVolume();
        };

        GridPoints grid; //convenient struct for storing gridpoints
        ScalarField rho; //density
        SymmetricTensorField Tv; //virial stress
        SymmetricTensorField Tk; //kinetic stress
        
    private:
        std::shared_ptr< ParticleSystem > m_sysdata; //!< particle system, equipped with neighbor list 
        std::shared_ptr< CoarseGrainFunction > m_cgfunc;    //!< coarse graining function
        //std::shared_ptr< PairPotential > ;         //!< pair potential
        double totTvxy; 
        double maxforce_rcut;
        double force_rcut;
        double cg_rcut;
        bool hasvelocities; 
};

//Compute Irving Kirkwood field
void IrvingKirkwood::compute()
{
    for(unsigned int n = 0; n < grid.totsize; n++)
    {
        for( auto p_i = euclidean_search(   m_sysdata->particles.get_query(), 
                                            vdouble3(grid.points[n][0],grid.points[n][1],grid.points[n][2]), 
                                            cg_rcut); p_i != false; ++p_i)
        {
            //Set up the grid position X and particle i position Ri in the cg function
            m_cgfunc->setX(grid.points[n]);
            Vector3d dr;
            dr << p_i.dx()[0], p_i.dx()[1],p_i.dx()[2];
            m_cgfunc->setRi(grid.points[n]+dr);
           
            //Safety check to see if dr < cg_rcut 
            if (dr.dot(dr) < cg_rcut*cg_rcut)
            {

                //Compute the density field
                double cg_val = m_cgfunc->getDeltaFunc();
                rho[n] += cg_val; 
                
                hasvelocities = false; 
                //If velocity of particles are included, we can compute the kinetic stress as well;
                if (hasvelocities)
                {
                    Tk.XX[n] += get<mass>(*p_i)*(get<velocity>(*p_i)[0])*(get<velocity>(*p_i)[0])*cg_val;
                    Tk.YY[n] += get<mass>(*p_i)*(get<velocity>(*p_i)[1])*(get<velocity>(*p_i)[1])*cg_val;
                    Tk.XY[n] += get<mass>(*p_i)*(get<velocity>(*p_i)[0])*(get<velocity>(*p_i)[1])*cg_val;
                    if (m_sysdata->simbox->getDim() > 2)
                    {
                        Tk.ZZ[n] += get<mass>(*p_i)*(get<velocity>(*p_i)[2])*(get<velocity>(*p_i)[2])*cg_val;
                        Tk.XZ[n] += get<mass>(*p_i)*(get<velocity>(*p_i)[0])*(get<velocity>(*p_i)[2])*cg_val;
                        Tk.YZ[n] += get<mass>(*p_i)*(get<velocity>(*p_i)[1])*(get<velocity>(*p_i)[2])*cg_val;
                    } 
                }

                //Next we loop through the j-th particles for the virial stress
                for( auto p_j = euclidean_search(   m_sysdata->particles.get_query(), 
                                                    get<position>(*p_i),maxforce_rcut); p_j != false; ++p_j)
                {
                    //Make sure the particle is unique
                    if (get<id>(*p_i) != get<id>(*p_j))
                    {
                        //set the distance between particle i and particle j
                        Vector3d rij;
                        rij << p_j.dx()[0], p_j.dx()[1], p_j.dx()[2];
                        m_cgfunc->setRij(rij);
                        m_sysdata->potential->setRij(rij);

                        //Don't forget to set diameters of the potential
                        m_sysdata->potential->setDiameters(get<diameter>(*p_i),get<diameter>(*p_j));
                        
                        //Make sure that it's within range.
                        if (rij.dot(rij) < 0.25*force_rcut*force_rcut*(get<diameter>(*p_i)+get<diameter>(*p_j)))
                        {
                            //Compute pair-force!
                            Vector3d F;
                            F = (m_sysdata->potential->getPairForce())*rij;//C.noalias() = A * B.transpose();
                            double bond_val = m_cgfunc->getBondFunc();
                            Tv.XX[n] += 0.5*F[0]*rij[0]*bond_val;
                            Tv.YY[n] += 0.5*F[1]*rij[1]*bond_val;
                            Tv.XY[n] += 0.5*F[0]*rij[1]*bond_val;
                            if (get<id>(*p_i) < get<id>(*p_j) && n == 0)
                            {
                                totTvxy += F[0]*rij[1];
                            }
                            
                            if (m_sysdata->simbox->getDim() > 2)
                            {
                                Tv.ZZ[n] += 0.5*F[2]*rij[2]*bond_val;
                                Tv.XZ[n] += 0.5*F[0]*rij[2]*bond_val;
                                Tv.YZ[n] += 0.5*F[1]*rij[2]*bond_val;
                            }
                        }
                    }
                }
            }
        } //end of for loop n
    }             
};

//Compute Irving Kirkwood field
void IrvingKirkwood::globalcompute()
{
    //#pragma omp parallel for reduction(+:totTvxy)
    for( auto p_i = m_sysdata->particles.begin(); p_i != m_sysdata->particles.end(); ++p_i)
    {
        //Next we loop through the j-th particles for the virial stress
        for( auto p_j = euclidean_search(   m_sysdata->particles.get_query(), 
                                            get<position>(*p_i),maxforce_rcut); p_j != false; ++p_j)
        {
            //Make sure the particle is unique
            if (get<id>(*p_i) < get<id>(*p_j))
            {
                //set the distance between particle i and particle j
                Vector3d rij;
                rij << p_j.dx()[0], p_j.dx()[1], p_j.dx()[2];
                m_sysdata->potential->setRij(rij);
                //Don't forget to set diameters of the potential
                m_sysdata->potential->setDiameters(get<diameter>(*p_i),get<diameter>(*p_j));
                
                //Compute pair-force!
                double forceval = m_sysdata->potential->getPairForce();
                Vector3d F =  (forceval)*rij;
                totTvxy += 0.5*(F[0]*rij[1]+F[1]*rij[0]);
            }
        }
    } //end of for loop n
};

void export_IrvingKirkwood(py::module& m)
{
    py::class_<IrvingKirkwood, std::shared_ptr<IrvingKirkwood> >(m,"IrvingKirkwood")
    .def(py::init< std::vector< Vector3d >, std::shared_ptr< ParticleSystem >, std::shared_ptr< CoarseGrainFunction > >())
    .def("compute", &IrvingKirkwood::compute)
    .def("globalcompute", &IrvingKirkwood::globalcompute)
    .def("getRho", &IrvingKirkwood::getRho)
    .def("getTv_xy", &IrvingKirkwood::getTv_xy)
    .def("getTvxy", &IrvingKirkwood::getTvxy)
    ;
};
