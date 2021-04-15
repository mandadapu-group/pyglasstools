import pyglasstools.utils as utils
import pyglasstools.io as io
import pyglasstools.potential as pot
import pyglasstools as pglass
import pyglasstools.elasticity as els
import pyglasstools.thermo as thermo
import numpy as np

start = 0
lj = pot.polydisperse10(v0=1.0,eps=0.1,rcut=1.4)
mysystem = utils.read_gsd(filename='inherent.gsd',checkpointfile="restart.log",ndim=3)

#Define nonaffine tensor calculator
myeigensolver = els.eigensolver(package='slepc-petsc')
names = ['nonaffinetensor_0000','nonaffinetensor_0011','nonaffinetensor_0001','nonaffinetensor_1111','nonaffinetensor_1101','nonaffinetensor_0101']
mylogger = els.logfile(filename="newnonaffine.log",names=names,solver=myeigensolver)
names = ['nconv']
for i in range(10):
    names.append('eigenvalue_{}'.format(i))
    names.append('eigenrelerror_{}'.format(i))
mylogger = els.logfile(filename="eigensummary.log",names=names,solver=myeigensolver)

#names_list = ['eigenvector_0','eigenvector_1','eigenvector_2','eigenvector_3']
#mylogger1 = els.fieldlogger(keyword="eigen", names=names_list,solver=myeigensolver)
mythermocalculator = thermo.calculator()
names = ['borntensor_0000','borntensor_0011','borntensor_0001','borntensor_1111','borntensor_1101','borntensor_0101']
mylogger1 = io.logfile(filename="newborntensor.log",names=names,solver=mythermocalculator)
names = ['virialstress_00','virialstress_11','virialstress_01']
mylogger2 = io.logfile(filename="inherent-virialstress.log",names=names,solver=mythermocalculator)

pglass.analyze(frame_list = range(45))#len(mysystem.traj)))
