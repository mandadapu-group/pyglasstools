# **PyGlassTools**

PyGlassTools is a Python module that compiles a couple of calculations for analyzing atomistic trajectories that are routinely done in our group, particularly supercooled liquids and glasses. Flexibility is important we don't deal not just with one type of pair potential energy function, but multiple kinds that represent various glass-forming liquids. Most of these calculations are typically Irving-Kirkwood calculations, normal mode analysis, and elasticity tensor calculations. 

The package is ready to use, but bear with us until more detailed instructions on installation and tutorials will come! 

## **Announcements**

(03/08/2021): Currently in the process of re-factoring as well as moving functions and class methods to separate implementation files to reduce (re)compilation time.

## **Contents** 

Files that come with this module:
 - CMakeLists.txt           : main CMake configuration file for the plugin
 - CompilerFlagsSetup.cmake : CMake file for setting up compiler flags
 - SetupPybind11.cmake      : CMake file for setting Python build using [pybind11](https://pybind11.readthedocs.io/en/stable/)
 - README.md                : This file
 - pyglasstools             : Directory containing C++ and Python source codes that make up the module
 - setup.py                 : Python script to install PyGlassTools as a module. *DO NOT RUN IT* on the source directory. 

## **Requirements**

For general use:
- Python >= 3.5
- CMake >= 2.8.10.1
- gcc and g++, capable of C++14 (tested with gcc-5, gcc-7, and gcc-9)
- [Boost](https://www.boost.org/) >= 1.50, including its [serialization library](https://www.boost.org/doc/libs/1_72_0/libs/serialization/doc/index.html) 
- [OpenMPI](https://www.open-mpi.org/) >=1.10

For irvingkirkwood module
- [GSL](https://www.gnu.org/software/gsl/) >= 2.0 (required for Irving-Kirkwood calculations)

For elasticity module
- [PETSc](https://www.mcs.anl.gov/petsc/) (with [ScaLAPACK](http://www.netlib.org/scalapack/), [MUMPS](http://mumps.enseeiht.fr/), PETSc's bundled [BLAS/LAPACK](https://bitbucket.org/petsc/pkg-fblaslapack/src/master/) installed "internally").
- [SLEPc](https://slepc.upv.es/)

## **Installation Instructions**

Before cloning the project, please install all the requirements! At the moment, you need to install all requirements to use the package (even if you are interested in only using one of our submodules). With the exception to PETSc/SLEPc, the installation of most of these requirements are fairly straightforward. 

### **Notes on SLEPc/PETSc Installation**

SLEPc/PETSc/ is a parallel linear algebra (PETSc) and eigensolver (SLEPc) work-horse. Its installation comes with multiple customizing options. Prior to installing SLEPc, you should have PETSc installed first!

For PETSc installation, we ask you to add the following additional flags to their configure file:
```console
--download-fblaslapack --download-scalapack --download-mumps
```
These are requirements are necessarry for the specific usage of the elasticity module in pyglasstools. (see **Example Scripts** section) 

### **Install Package**

First, git clone the project:
```console
$ git --recursive clone https://github.com/mandadapu-group/pyglasstools
```
Don't forget the recursive flag since this project utilizes git submodules. 

Next, build the project
```console
$ cd pyglasstools
$ mkdir build
$ cd build
$ cmake ../ 
$ make -j4
```

Afterwards, install with pip
```console
$ pip install .
```

## **Example Scripts**

We list several examples

## **To-Do List**

Features that are really should be there, but don't have enough time at the moment:
1. Implementation of arbitrarily-defined user pair potentials.
2. Implementation of arbitrarily-defined user observables (both macroscopic and IK fields).
3. Alternative package for eigendecomposition analysis (besides SLEPc, which may be too hard to use!).

Features that are relatively desirable, but not in high demand at the moment:
1. Implement parallelized excitation analysis (in the context of Dynamical Facilitation Theory) 
2. (*Optional*) Implement stress and density auto-correlation functions
3. (*Optional*) Interface for 'quick' plotting and movie generation. 
