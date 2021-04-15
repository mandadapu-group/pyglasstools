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
