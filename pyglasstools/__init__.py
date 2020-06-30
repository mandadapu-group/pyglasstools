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
from tqdm import tqdm


_default_excepthook = sys.excepthook;
## \internal
# \brief Override pythons except hook to abort MPI runs
def _pyglasstools_sys_excepthook(type, value, traceback):
    _default_excepthook(type, value, traceback);
    sys.stderr.flush();
    comm.abort();

sys.excepthook = _pyglasstools_sys_excepthook

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
    progressbar = tqdm(total=len(frame_list),file=sys.stdout)
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
                logger.save(frame_num, savemode);
        if rank == 0:
            progressbar.update(1)
            print("")
        comm.barrier()

def reset():
    global loggers_list, solvers_list, savemode
    loggers_list.clear()
    solvers_list.clear()
    savemode = "cartesian"
