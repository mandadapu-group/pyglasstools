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
