

# **PyGlassTools**

PyGlassTools is a Python module that compiles a couple of calculations for analyzing atomistic trajectories that are routinely done in our group, particularly ones for supercooled liquids and glasses. Flexibility is important since we don't deal with just with one type of pair potential energy function, but multiple kinds that represent various glass-forming liquids. Most of these calculations are typically Irving-Kirkwood calculations, normal mode analysis, and elasticity tensor calculations. 

The package is primarily used for all the analysis used in this paper:
- M. R. Hasyim, K. K. Mandadapu, "A Theory of Localized Excitations in Supercooled Liquids" https://arxiv.org/abs/2103.03015 (2021)

The package is ready to use for 2D systems, but bear with us until more detailed instructions on installation and tutorials will come! 

If you do have any issues, bring it up on GitHub or contact me via e-mail (muhammad_hasyim@berkeley.edu)

## **Announcements**

(05/02/2021): Still in the process of re-factoring as well as moving functions and class methods to separate implementation files to reduce (re)compilation time. Elasticity module now works in 3D. 

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
- [GSD](https://gsd.readthedocs.io/en/stable/index.html)

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
These requirements are necessarry for the specific usage of the elasticity module in pyglasstools. (see **Example Scripts** section) 

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

## **Introduction to PyGlassTools**


### **Please Read Before Running**

The package only reads .gsd trajectory files. This trajectory file is the preferred file format when using the MD/MC simulation package [HOOMD-Blue](https://github.com/glotzerlab/hoomd-blue). Any trajectory file you have must be converted to .gsd files prior to using this package. See the [GSD home page](https://gsd.readthedocs.io/en/stable/index.html) for installation and examples on how to use the package. 

Unfortunately, we also have yet to implement pair potentials that can be arbitrarily defined by the user. At the moment, there are pre-compiled pair potentials available to use. Nearly all pair potentials are geared solely for the analysis of atomistic glass formers with continuous poly-dispersity. The model is defined by the following pair potential

![equation](https://latex.codecogs.com/gif.latex?%5Cphi%28r/%5Csigma_%7B%5Calpha%20%5Cbeta%7D%29%20%3D%20%5Cbegin%7Bcases%7D%20v_0%20%5Cleft%5B%5Cleft%28%5Cdfrac%7B%5Csigma_%7B%5Calpha%20%5Cbeta%7D%7D%7Br%7D%5Cright%29%5Em-%5Cleft%28%5Cdfrac%7B%5Csigma_%7B%5Calpha%20%5Cbeta%7D%7D%7Br%7D%5Cright%29%5En%5Cright%5D&plus;%5Csum_%7Bk%3D0%7D%5Eq%20c_k%20%5Cleft%28%5Cfrac%7Br%5E%7B%5Calpha%20%5Cbeta%7D%7D%7B%5Csigma_%7B%5Calpha%20%5Cbeta%7D%7D%20%5Cright%20%29%5E%7B2k%7D%26%20r/%5Csigma_%7B%5Calpha%20%5Cbeta%7D%20%5Cleq%20%5Ctilde%7Br%7D_c%20%5C%5C%200%20%26%20%5Ctext%7Botherwise%7D%20%5Cend%7Bcases%7D)

![equation](https://latex.codecogs.com/gif.latex?%5Csigma_%7B%5Calpha%20%5Cbeta%7D%20%3D%20%5Cfrac%7B1%7D%7B2%7D%5Cleft%28%5Csigma_%5Calpha%20&plus;%5Csigma_%5Cbeta%5Cright%29%281-%5Cvarepsilon%7C%5Csigma_%5Calpha%20-%20%5Csigma_%5Cbeta%7C%29)

The first term in the first equation is the standard repulsive and attractive interaction. The second term is an even polynomial ensuring smoothness up to q-th order at the cut off radius. 

The models available to use and currently implemented are as follows:

|   Model Name      |   q       |   m       |   n       | 
|   :--------       |   :--:    |   :--:    |   :--:    |
|   polydisperse12  |   2       |   12      |   0       |
|   polydisperse18  |   2       |   18      |   0       |
|   polydisperselj  |   2       |   12      |   6       |
|   polydisperse10  |   3       |   10      |   0       |
|   polydisperse106 |   2       |   10      |   6       |

You will see in pyglasstools/potential/\_\_init\_\_.py file that there are other pair potentials, but I haven't thoroughly tested them! So be aware. 

Lastly, some module only works in 2D while others work both in 2D and 3D. For instance, we have not extended the irvingkirkwood module to 3D.

All scripts shown here are part of the example_scripts of this project! These test scripts are applied to the polydisperse12 trajectory file. This trajectory consists of 201 energy-minimizing configurations of N=32^2 particles.  

### **pyglasstools.thermo** 

This module computes various macroscopic observables of the system given a pair potential of the system. Most MD packages already have that feature, but there are lesser known observables (like the Born part of the elasticity tensor) which is not typically implemented. 

Example script is located in pyglasstools/example_scripts/thermo 

script name: *test-thermo-module.py*
```python
import pyglasstools as pglass
import pyglasstools.utils as utils
import pyglasstools.io as io
import pyglasstools.potential as pot
import pyglasstools.thermo as thermo
#import numpy as np

#Construct the pair potential
poly12 = pot.polydisperse12(v0=1.0,eps=0.2,rcut=1.25)

#Read the system from a GSD file. Also construct checkpoint file restart.log, in case the script stops unexpectedly
mysystem = utils.read_gsd(filename='../polydisperse12.gsd',checkpointfile="restart.log") 

#Construct the a calculator object, which computes all observables that we want to compute
MyCalculator = thermo.calculator()

#For this example, we will compute two different observables. The Born part of the elasticity tensor and the virial stress. 
#First we do the Born elasticity tensor. Index runs from 0,1,2 indicating x,y, and z components respectively.
names = ['borntensor_0000','borntensor_0011','borntensor_0001','borntensor_1111','borntensor_1101','borntensor_0101']
mylogger = io.logfile(filename="borntensor.log",names=names,solver=MyCalculator)

#Second, we do the virialstress
names = ['virialstress_00','virialstress_11','virialstress_01']
mylogger1 = io.logfile(filename="virialstress.log",names=names,solver=MyCalculator)

#Now, we run the analysis!
pglass.analyze(frame_list = range(200))
```

Once, we run this script, we can check one interesting property of the system which is isotropy. In an isotropic system, the shear modulus can be computed from various components of the elasticity tensor. We can quickly check to see if this is true!

script name: *check-isotropy.py*
```python
import numpy as np

#Since the system is isotropic, shear moduli can be computed from various components of the tensor. 
borntensor = np.loadtxt('borntensor.log',skiprows=1)

#First component is the xyxy component
shearmod1 = np.mean(borntensor[:,6])

#Second component comes from the xxxx and xxyyy component
shearmod2 = 0.5*np.mean(borntensor[:,1]-borntensor[:,2])

#Check if they're the same up some significant figures
print(shearmod1,shearmod2)
```

The output should be
```console
$ 14.89655895 14.8769594875
```
where we see that the Born shear modulus computed in two ways agree by 3 significant figures. More samples and larger system sizes will make this result even more accurate!

### **pyglasstools.irvingkirkwood**

This module computes Irving-Kirkwood fields of various observables (the same ones you see in __pyglasstools.thermo__). For example, the virial stress:

![equation](https://latex.codecogs.com/gif.latex?%5Cmathbf%7BT%7D%28%5Cmathbf%7Bx%7D%29%20%3D%20%5Csum_%7B%5Calpha%20%5Cbeta%7D%5Cmathbf%7BF%7D_%7B%5Calpha%20%5Cbeta%7D%20%5Cotimes%20%5Cmathbf%7Br%7D_%7B%5Calpha%20%5Cbeta%7D%20%5Cint_0%5E1%20%5Ctext%7Bd%7Ds%20%5C%20%5CDelta%28%5Cmathbf%7Bx%7D-%5Cmathbf%7Br%7D_%7B%5Calpha%7D&plus;s%5Cmathbf%7Br%7D_%7B%5Calpha%20%5Cbeta%7D%29)

The module can compute these stress fields with using a selection of coarse graining functions and quadrature methods (to compute the integral).

script name: *test-irvingkirkwood-module.py*
```python
import pyglasstools as pglass  
import pyglasstools.io as io  
import pyglasstools.utils as utils  
import pyglasstools.potential as pot  
import pyglasstools.irvingkirkwood as ik  
import pyglasstools.irvingkirkwood.cgfunc as cgfunc  
import numpy as np  
import timeit  
  
#Construct the pair potential  
poly12 = pot.polydisperse12(v0=1.0,eps=0.2,rcut=1.25)  
  
#Construct the coarse-graining function  
mycgfunc = cgfunc.octic(rcut = 2.0, mode="fixedpoint",order=5)  
  
#Read the system from a GSD file. Also construct checkpoint file restart.log, in case the script stops unexpectedly  
mysystem = utils.read_gsd(filename='../polydisperse12.gsd',checkpointfile="restart.log")  
  
#Construct the calculator object, which computes the coarse grined fields we want to compute. dx sets the spacing of the lattice grid we compute the stress with.  
MyCalculator = ik.ikcalculator(cgfunc=mycgfunc,dx=0.25)  
  
#Let us compute the virial stresses  
names = ['virialstress_00','virialstress_11','virialstress_01']  
mylogger = io.fieldlogger(keyword="test-dump",names=names,solver=MyCalculator)  
pglass.analyze(frame_list = range(180,200))
```

The module can be run normally
```console
$ python test-irvingkirkwood-module.py
```
 or in parallel
```console
$ mpirun -n 5 python test-irvingkirkwood-module.py
```
Certain situations calls for parallel runs (large system sizes). In our case, this is not needed. 

Running this module requires the help of a submodule called __pyglasstools.irvingkirkwood.cgfunc__ this module contains various coarse graining functions to use for the computation. 

|	Coarse Graining Function	|	Equation 	|
|	:--	|	:--: |
|	octic | ![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7B15%7D%7B8%20%5Cpi%20%28r_%5Cmathrm%7Bcut%7D%29%5E2%7D%5Cleft%281-2%20%5Cleft%28%5Cfrac%7Br%7D%7Br_%5Cmathrm%7Bcut%7D%7D%20%5Cright%20%29%5E4%20&plus;%20%5Cleft%28%5Cfrac%7Br%7D%7Br_%5Cmathrm%7Bcut%7D%7D%20%5Cright%20%29%5E8%20%5Cright%29) |
|rect	| ![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7B1%7D%7B%5Cpi%28r_%5Cmathrm%7Bcut%7D%29%5E2%7D)|
|mollifier | ![equation](https://latex.codecogs.com/gif.latex?%5Cfrac%7BC_1%7D%7B2%20%5Cpi%20%28r_%5Cmathrm%7Bcut%7D%29%5E2%7D%5Cexp%5Cleft%28-%5Cfrac%7B1%7D%7B1-%20%5Cleft%28%5Cfrac%7Br%7D%7Br_%5Cmathrm%7Bcut%7D%7D%20%5Cright%20%29%5E2%20%7D%5Cright%29%3B%20%5Cquad%20C_1%20%3D13.468420987%5Cldots) |

For all three coarse graining functions, we can choose two different integration methods/modes. This is chosen from the 'mode' argument when we first construct the coarse-graining function in the script. Implementation of both methods are from GSL. 

|  Integration Mode | Description | Parameters|
|	:----		| :--- |  :----|
| 'fixedpoint' |  Gauss-Legendre integration | _order_:  the order of quadrature |
| 'adaptive' | QAG adaptive integration | _relerr_: relative error, _abserr_: absolute error, _maxiter_: maximum iterations |

For more details, see GSL's documentation [here](https://www.gnu.org/software/gsl/doc/html/integration.html). Note that 'fixedpoint' requires no iteration when computing the integral, and thus it is the fastest method of computing the integral. 

### **pyglasstools.elasticity**

Despite its name, this module is targeted specifically for the computation of the non-affine part of the elasticity tensor (the Born part is already handled by the __pyglasstools.thermo__  module). Because its computation requires the computation of all normal modes of the system, it also gives us information about vibrational density of states in the process.  

The work-horse of this module is SLEPc, the parallel eigensolver library. In fact, one may think of this module as a small wrapper to SLEPc. Thus, familiarity with how SLEPc is used in a simpler code will be a huge help (see the manual [here](https://slepc.upv.es/documentation/slepc.pdf]) )

In the example script below, we demonstrate how this module can be used in conjunction with __pyglasstools.thermo__ to compute both the Born and non-affine part of the elasticity tensor!

script name: _test-elasticity-module.py_
```python
import pyglasstools.utils as utils  
import pyglasstools.io as io  
import pyglasstools.potential as pot  
import pyglasstools as pglass  
import pyglasstools.elasticity as els  
import pyglasstools.thermo as thermo  
import numpy as np  
  
poly12 = pot.polydisperse12(v0=1.0,eps=0.2,rcut=1.25)  
mysystem = utils.read_gsd(filename='../polydisperse12.gsd',checkpointfile="restart.log",ndim=2)  
  
#Construct the eigensolver and choose which package to use as backend for sigensolver  
myeigensolver = els.eigensolver(package='slepc-mumps')  
  
#We want to compute the non-affine part of the elasticity tensor  
names = ['nonaffinetensor_0000','nonaffinetensor_0011','nonaffinetensor_0001','nonaffinetensor_1111','nonaffinetensor_1101','nonaffinetensor_0101']  
mylogger = els.logfile(filename="newnonaffine.log",names=names,solver=myeigensolver)  
  
#In addition, we want to store the total number of converged eigenmodes 'nconv',  
#the first 10 eigenvalues found during the computation and their relative error  
names = ['nconv']  
for i in range(10):  
names.append('eigenvalue_{}'.format(i))  
names.append('eigenrelerror_{}'.format(i))  
mylogger = els.logfile(filename="eigensummary.log",names=names,solver=myeigensolver)  
  
#We can also store the first four eigenmodes found by this package  
names_list = ['eigenvector_0','eigenvector_1','eigenvector_2','eigenvector_3']  
mylogger1 = els.fieldlogger(keyword="eigen", names=names_list,solver=myeigensolver)  
  
#Next, we construct the calculator used to compute global observables like the Born elasticity tensor and virial stress.  
mythermocalculator = thermo.calculator()  
  
names = ['borntensor_0000','borntensor_0011','borntensor_0001','borntensor_1111','borntensor_1101','borntensor_0101']  
mylogger1 = io.logfile(filename="newborntensor.log",names=names,solver=mythermocalculator)  
  
names = ['virialstress_00','virialstress_11','virialstress_01']  
mylogger2 = io.logfile(filename="inherent-virialstress.log",names=names,solver=mythermocalculator)  
  
#Let's analyze the first two frames for now!  
pglass.analyze(frame_list = range(2))
```

To run the script, additional command lines are needed. Here's an example:

```console
$ mpirun -np 8 python test-elasticity-module.py -lowerbound_tol -1.0 -upperbound_tol 1e-6 -pinv_tol 1e-6 -eps_tol 1e-7 -eps_max_it  
10000 -notice_level 3
```
The following is explanation for every command-line option:
| Command Line Option | Description |
| :---: | :-- |
| -lowerbound_tol | Tolerance for the lowest eigenvalue we want to compute |
| -upperbound_tol | Small tolerance above the highest eigenvalue of the system to ensure no false detection |
| -pinv_tol | Tolerance for the smallest magnitude of detected eigenvalues, before they're considered numerically zero |
| -eps_tol | Tolerance for convergence of each eigenmode |
| -eps_max_it | Maximum number iteration for the eigensolver |
| -notice_level | The level of warnings and print messages during the computation |

To compute the non-affine part of the elasticity tensor, we want to make sure that all eigenmodes are computed. This is why we set the lowest eigenvalue to be a negative value (since the .gsd file included by the example_scripts contains configuration of particles whose Hessian is semi-positive definite). 

Among these command-line options, only -eps_tol and -eps_max_it are are used by SLEPc back-end package. The rest are input options taken in independently by pyglasstools for configuring the eigensolver. Thus, when running the script, you will encounter the following warning message:

```console
WARNING! There are options you set that were not used!  
WARNING! could be spelling mistake, etc!  
There are 4 unused database options. They are:  
Option left: name:-lowerbound_tol value: -1.0  
Option left: name:-notice_level value: 3  
Option left: name:-pinv_tol value: 1e-6  
Option left: name:-upperbound_tol value: 1e-6
``` 
This warning is an automatic message from SLEPc back-end warning messages saying that it did not use these options, i.e., they are not part of SLEPc's database of command-line options. However, in reality, pyglasstools does use these options. We will fix this issue in the near future so it won't mislead users!

Note that the eigensolver has a 'package' argument, indicating which back-end we use for solver. SLEPc is unique in a sense that it can rely on various linear algebra packages that PETSc has available for you. Currently, we let user to choose two packages:

| Eigensolver Package Name | Description|
|  :---: | :-- |
| 'slepc-mumps' | Uses the parallel linear algebra package as [MUMPS](http://mumps.enseeiht.fr/) and [ScaLAPACK](http://www.netlib.org/scalapack/) |
| 'slepc-petsc' | Uses PETSc's internal linear algebra solver |

If 'slepc-petsc' is chosen, then the an additional command-line must be given called -eps_krylovschur_partitions. The argument of this option must match the number of MPI processes used:
```console
$ mpirun -np 8 python test-elasticity-module.py -lowerbound_tol -1.0 -upperbound_tol 1e-6 -pinv_tol 1e-6 -eps_tol 1e-7 -eps_max_it  
10000 -notice_level 3 -eps_krylovschur_partitions 8
```
The reasoning behind it is a bit technical, and it deals with how PETSc's linear algebra solver is parallelized for the purpose of eigensolving. Details of this can be found in SLEPc manual [here](https://slepc.upv.es/documentation/slepc.pdf]). We recommend using the 'slepc-mumps' because it is more robust against a common problem in PETSc known as [the zero pivot problem](https://www.mcs.anl.gov/petsc/documentation/faq.html#zeropivot).

## **To-Do List**

Features that are really should be there, but don't have enough time at the moment:
1. Arbitrary input files (not just .gsd files, but also the more common .xyz files) 
2. Implementation of arbitrarily-defined user pair potentials.
3. Implementation of arbitrarily-defined user observables (both macroscopic and IK fields).
4. Alternative package for eigendecomposition analysis (besides SLEPc, which may be too hard to use!).
5. Doxygen for C++ classes documentation.

Features that are relatively desirable, but not in high demand at the moment:
1. Switching the custom MPI interface to [MPL](https://github.com/rabauke/mpl).
2. Implement parallelized [excitation analysis](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.1.021013) (in the context of Dynamical Facilitation Theory). 
3. Implement stress and density auto-correlation function calculations.
4. An interface for 'quick' plotting and movie generation. 
