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

## **To-Do List**
1. Local calculator and stress field calculation
2. Elasticity tensor calculation
3. Eigenvalue solvers for many tasks:
   * Computing pseudo-inverse of the Hessian matrix
   * Computing the soft-modes of the glass
4. (*Optional*) Implement stress and density auto-correlatoin functions
5. (*Optional*) Interface for 'quick' plotting and movie generation. 

Note that for the eigensolver, this will require choosing an external library. Multiple options to use different external library would be highly desirable. Example: [Spectra](https://spectralib.org/), [SLEPc](https://slepc.upv.es/), and [FEAST](http://www.ecs.umass.edu/~polizzi/feast/). 

For the pseudoinverse problem, libraries specializing on (truncated) SVD would be extremely useful.
