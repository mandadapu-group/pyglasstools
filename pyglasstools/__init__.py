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

#Initialize a single communicator during module call
comm = _pyglasstools.Communicator()
rank = comm.getRank()
size = comm.getNRanks()

#Save mode
savemode = "unordered"
def set_savemode(name):
    savemode = name

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
globalsysdata = None
globalpotential = None

def set_sysdata(sysdata):
    global globalsysdata
    globalsysdata = sysdata

def update_sysdata(frame_num):
    global globalsysdata
    globalsysdata.update(frame_num)

def get_sysdata():
    global globalsysdata
    return globalsysdata

def set_potential(potential):
    global globalpotential
    globalpotential = potential

def get_potential():
    global globalpotential
    return globalpotential


import atexit
#This is an important function for anyone who wants to tinker with how the analysis is run and also for "cleanly" exiting the module
def reset():
    global loggers_list, solvers_list
    loggers_list.clear()
    solvers_list.clear()
atexit.register(reset)

import time
def wait_until(timeout, period=0.25):
    global ioflag
    mustend = time.time() + timeout
    while time.time() < mustend:
        if ioflag: 
            return True
        time.sleep(period)
    return False

import signal
import logging
import logging.handlers
log = logging.getLogger()

# class based on: http://stackoverflow.com/a/21919644/487556
class DelayedInterrupt(object):
    def __enter__(self):
        self.signal_received = False
        self.old_handler = {}
        self.old_handler[signal.SIGINT] = signal.signal(signal.SIGINT, self.handler)
        self.old_handler[signal.SIGTERM] = signal.signal(signal.SIGTERM, self.handler)

    def handler(self, sig, frame):
        self.signal_received = (sig, frame)
        print(f'Signal {sig} received by Process [{rank}]. Delaying KeyboardInterrupt.')
    def __exit__(self, type, value, traceback):
        for sig in [signal.SIGINT, signal.SIGTERM]: 
            signal.signal(sig, self.old_handler[sig])
            if self.signal_received:
                self.old_handler[sig](*self.signal_received)

from pyglasstools import utils;
from random import randint, uniform

def analyze(frame_list):
    global globalsysdata
    totallength = len(frame_list)
    frame_list = globalsysdata.setup_checkpoint(frame_list)
    if not frame_list:
        if rank == 0:
            print("No more frames to analyze. Check your logfiles/output")
    else:
        if rank == 0:
            progressbar = tqdm(total=totallength,file=sys.stdout,leave=False,initial=frame_list[0])
            print("")
        
        for frame_num in frame_list:
            for solver in solvers_list:
                solver.update(frame_num);
                solver.run();
            with DelayedInterrupt():
                for logger in loggers_list:
                    logger.save(frame_num);
                #It's sufficient to go use at least one solver to do the checkpointing
                if rank == 0 and globalsysdata.checkpointfile is not None:
                    with open(globalsysdata.checkpointfile,"w") as f:
                        f.write("Frame {}".format(frame_num+1))
                        f.close()
                if rank == 0:
                    progressbar.update(1)
                    print("")
                #Barrier for consistent update
                if size > 1:
                    comm.barrier()
                pass 
        if rank == 0:
            progressbar.close()
