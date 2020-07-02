""" PyGlassTools python API
"""
import sys;
import ctypes;
import os;
from tqdm import tqdm
import numpy as np

#In order for our MPI wrapper to work , we need to set dlopen flag to need to RTLD_GLOBAL 
flags = sys.getdlopenflags();
sys.setdlopenflags(flags | ctypes.RTLD_GLOBAL);

from pyglasstools import _pyglasstools;
from pyglasstools import utils;

#Initialize a single communicator during module call
comm = _pyglasstools.Communicator()
rank = comm.getRank()
size = comm.getSizeGlobal()

## \internal
_default_excepthook = sys.excepthook;

# \brief Override pythons except hook to abort MPI runs
def _pyglasstools_sys_excepthook(type, value, traceback):
    _default_excepthook(type, value, traceback);
    sys.stdout.flush();
    sys.stderr.flush();
    comm.abort(1);

sys.excepthook = _pyglasstools_sys_excepthook


#The module should save more than one logger, which analyzes various observables
loggers_list = [];
solvers_list = [];

import atexit
#This is an important function for anyone who wants to tinker with how the analysis is run and also for "cleanly" exiting the module
def reset():
    global loggers_list, solvers_list
    loggers_list.clear()
    solvers_list.clear()
atexit.register(reset)


def analyze(frame_list,mode="normal"):
    if rank == 0:
        progressbar = tqdm(total=len(frame_list),file=sys.stdout,leave=False,initial=frame_list[0])
        print("")
    for frame_num in frame_list:
        for solver in solvers_list:
            solver.update(frame_num);
            if (mode == "normal"):
                solver.run();
            else:
                solver.run(mode);
        for logger in loggers_list:
            logger.save(frame_num);
        if rank == 0:
            progressbar.update(1)
            print("")
    if rank == 0:
        progressbar.close()
