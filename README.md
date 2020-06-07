# **PyGlassTools**

PyGlassTools is a Python module that compiles a couple of calculations that I routinely do when analyzing atomistic trajectories, particularly supercooled liquids and glasses. Flexibility is important as I deal not just with one type of potential energy function, but multiple kinds that represent various glass-forming liquids. Most of these calculations are typically Irving-Kirkwood calculations, normal mode analysis, and elasticity tensor calculations. 

This is a (highly) experimental version for development purposes. **Do not git pull/clone** unless you are confident with what you're doing!

## **Announcements**

(06/06/2020): Most of the restructuring is finished. Added force dipole calculations and a way for nonaffine and irvingkirkwood modules to communicate. Will do more test calculations soon.

## **Contents** 

Files that come with this module:
 - CMakeLists.txt           : main CMake configuration file for the plugin
 - CompilerFlagsSetup.cmake : CMake file for setting up compiler flags
 - SetupPybind11.cmake      : CMake file for setting Python build using [pybind11](https://pybind11.readthedocs.io/en/stable/)
 - README.md                : This file
 - pyglasstools             : Directory containing C++ and Python source codes that make up the module
 - setup.py                 : Python script to install PyGlassTools as a module. *DO NOT RUN IT*. Special instructions will come.

## **Requirements**

For general use:
- Python >= 3.5
- CMake >= 2.8.10.1
- gcc and g++, capable of C++14 (tested with gcc-5, gcc-7, and gcc-9)
- [Boost](https://www.boost.org/) >= 1.50, including its [serialization library](https://www.boost.org/doc/libs/1_72_0/libs/serialization/doc/index.html) 
- [OpenMPI](https://www.open-mpi.org/) >=1.10

For irvingkirkwood module
- [GSL](https://www.gnu.org/software/gsl/) >= 2.0 (required for Irving-Kirkwood calculations)

For nonaffine module
- [PETSc](https://www.mcs.anl.gov/petsc/) (with [ScaLAPACK](http://www.netlib.org/scalapack/), [MUMPS](http://mumps.enseeiht.fr/), PETSc's bundled [BLAS/LAPACK](https://bitbucket.org/petsc/pkg-fblaslapack/src/master/) installed "internally").
- [SLEPc](https://slepc.upv.es/)

## **Installation Instructions**
Coming soon . . .


## **To-Do List**
1. (*Optional*) Implement stress and density auto-correlation functions
2. (*Optional*) Interface for 'quick' plotting and movie generation. 

Current eigensolver uses [Spectra](https://spectralib.org/) and.or [SLEPc](https://slepc.upv.es/)
