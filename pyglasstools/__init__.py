""" PyGlassTools python API
"""
import sys;
import ctypes;
import os;

#In order for our MPI wrapper to work , we need to set dlopen flag to need to RTLD_GLOBAL 
flags = sys.getdlopenflags();
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL);

from pyglasstools import _pyglasstools;
from pyglasstools import utils;
from os import path
import numpy as np

#Initialize a single communicator during module call
comm = _pyglasstools.Communicator()
rank = comm.getRank()
size = comm.getSizeGlobal()

#The module should save more than one logger, which analyzes various observables
loggers_list = [];
solvers_list = [];
savemode = "cartesian"
def set_savemode(inmode):
    global savemode
    savemode = inmode

def analyze(frame_list,mode="normal"):
    for frame_num in frame_list:
        for solver in solvers_list:
            solver.update(frame_num);
            if (mode == "normal"):
                solver.run();
            else:
                solver.run(mode);
        comm.barrier()
        for logger in loggers_list:
            if savemode == "polar":
                logger.save(frame_num, savemode);
            else:
                logger.save(frame_num);

def reset():
    global loggers_list, solvers_list, savemode
    loggers_list.clear()
    solvers_list.clear()
    savemode = "cartesian"
