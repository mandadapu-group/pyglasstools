#ifndef __IRVING_KIRKWOOD_H__
#define __IRVING_KIRKWOOD_H__
#include "Observables.h"

class PYBIND11_EXPORT VirialStress : public Observable
{
    public:
        VirialStress(int _dim) : Observable("Virial Stress", "TENSOR", false, true, _dim){};
        ~VirialStress(){};
        
        void accumulate(    AboriaParticles::value_type particle_i, 
                            AboriaParticles::value_type particle_j, 
                            const std::shared_ptr<PairPotential>& potential)
            {
                
                //Compute pair-force!
                double forceval = potential->getPairForce();
                Eigen::Vector3d rij = potential->getRij();
                Eigen::Vector3d F =  (forceval)*(rij);
                
                //Compute virial stress
                val(0,0) += F[0]*(rij[0]);
                val(1,1) += F[1]*(rij[1]);
                val(0,1) += F[0]*(rij[1]);
                val(1,0) += F[1]*(rij[0]);
                
                if (dim > 2)
                {
                    val(2,2) += F[2]*(rij[2]);
                    val(0,2) += F[0]*(rij[2]);
                    val(2,0) += F[2]*(rij[0]);
                    val(1,2) += F[1]*(rij[2]);
                    val(2,1) += F[2]*(rij[1]);
                }
            }
};

class PYBIND11_EXPORT KineticStress : public Observable
{
    public:
        KineticStress(int _dim) : Observable("Kinetic Stress", "TENSOR", true, false, _dim){};
        ~KineticStress(){};
        
        void accumulate(AboriaParticles::value_type particle_i) 
            {
                val(0,0) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[0])*(abr::get<velocity>(particle_i)[0]);
                val(1,1) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[1])*(abr::get<velocity>(particle_i)[1]);
                val(0,1) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[0])*(abr::get<velocity>(particle_i)[1]);
                val(1,0) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[1])*(abr::get<velocity>(particle_i)[0]);
                if (dim > 2)
                {
                    val(2,2) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[2])*(abr::get<velocity>(particle_i)[2]);
                    val(0,2) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[0])*(abr::get<velocity>(particle_i)[2]);
                    val(2,0) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[2])*(abr::get<velocity>(particle_i)[0]);
                    val(1,2) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[1])*(abr::get<velocity>(particle_i)[2]);
                    val(2,1) += abr::get<mass>(particle_i)*(abr::get<velocity>(particle_i)[2])*(abr::get<velocity>(particle_i)[1]);
                }
            }
};

void export_VirialStress(py::module& m)
{
    py::class_<VirialStress, Observable, std::shared_ptr<VirialStress> >(m,"VirialStress")
    .def(py::init< int >()) 
    ;
};
void export_KineticStress(py::module& m)
{
    py::class_<KineticStress, Observable, std::shared_ptr<KineticStress> >(m,"KineticStress")
    .def(py::init< int >()) 
    ;
};

#endif
