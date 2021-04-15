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
