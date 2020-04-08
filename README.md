# **PyGlassTools**

PyGlassTools is a Python module that compiles a couple of calculations that I routinely do when analyzing atomistic trajectories, particularly supercooled liquids and glasses. Flexibility is important as I deal not just with one type of potential energy function, but multiple kinds that represent various glass-forming liquids.

This is a (highly) experimental version for development purposes. **Do not git pull/clone** unless you are confident with what you're doing! Current development in place is to construct a C++ backend that is also highly efficient

## **Contents** 

Files that come with this module:
 - CMakeLists.txt           : main CMake configuration file for the plugin
 - CompilerFlagsSetup.cmake : CMake file for setting up compiler flags
 - SetupPybind11.cmake      : CMake file for setting Python build using [pybind11](https://pybind11.readthedocs.io/en/stable/)
 - README.md                : This file
 - pyglasstools             : Directory containing C++ and Python source codes that make up the module

## **Requirements**

- Python >= 3.5
- NumPy >= 1.7
- CMake >= 2.8.10.1
- gcc and g++, capable of C++14 (tested with gcc-5 and gcc-9.2)
- [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [Boost](https://www.boost.org/)


## **Change Log**

(03/23/2020): Implemented simulation box class.

(03/24/2020): Implemented pair potential and coarse graining function classes.

(03/25/2020): Implemented system data class with neighbor list functionality using Aboria.

(03/27/2020): Implemented Irving Kirkwood class to compute global shear stres.

(03/28/2020): Changed module design and structure and included Observable class for general computations. 

(04/01/2020): Added bond stiffness calculation and local Irving Kirkwood calculation.

(04/02/2020): Improved speed on local Irving Kirkwood calculation and added small OpenMP features.

(04/03/2020): Added Born elasticity tensor calculation.

(04/06/2020): Added Hessian class and eigendecomposition analysis, using Spectra and SLEPc (through slepc4py).

## **To-Do List**
1. Computing pseudo-inverse of the Hessian matrix by eigendecomposition.
2. (*Optional*) Implement stress and density auto-correlatoin functions
3. (*Optional*) Interface for 'quick' plotting and movie generation. 

Current eigensolver uses [Spectra](https://spectralib.org/) and or [SLEPc](https://slepc.upv.es/), through its own Python binding, slepc4py. 
