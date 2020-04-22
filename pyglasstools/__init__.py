""" PyGlassTools python API
"""
from pyglasstools import _pyglasstools;
from pyglasstools import utils;
from pyglasstools import calculator;
from pyglasstools import observables
from os import path

#The module should save more than one logger, which analyzes various observables
loggers = [];

#will automatically create the required calculators . . .
class logger(object):
    global loggers 

    def __init__(self, filename=None, names = None, sysdata = None, pair = None):
        
        #Save filename
        self.filename = filename
        
        #First, parse the list of names based on what type of obsercables they are
        self.__glob_obs_names = [s for s in names if "g_" in s ];
        self.global_obs = [] 
        #self.__nonaffine_obs_names = [s for s in names if "na_" in s ];
        #self.nonaffine_obs = [] 
        
        #Initialize system data and pair potential of the system
        self.sysdata = sysdata;
        self.pair = pair;

        #Initialize the "calculators"
        self.ikglobal = None;
        
        if not (not self.__glob_obs_names):
            #Construct the global calculator
            self.ikglobal = _pyglasstools.GlobalCalculator(self.sysdata._getParticleSystem(),self.pair._getPairPotential())
            #Then, add observables
            self.global_obs = observables.initialize_global(self.__glob_obs_names)
            for name in self.global_obs:
                self.ikglobal.addObservable(self.global_obs[name]._getObservable())
        
        #Once initialization is done, the logger adds itself to the global list of available loggers 
        loggers.append(self)

        #Initialize file to save:
        self.file = open(self.filename,"w+")

        #Create column headers
        self.file.write("Frame ")
        for name in self.__glob_obs_names:
            self.file.write("{} ".format(name))
        self.file.write("\n")

    def run(self): 
        if not (not self.__glob_obs_names):
            self.ikglobal.compute()
    
    def save(self,frame_num):
        self.file.write("{:d} ".format(frame_num))
        for name in self.__glob_obs_names:
            val = 0
            if "g_virialstress" in name:
                val = self.global_obs['g_virialstress'].getVal()[int(name[-2]),int(name[-1])]
            elif "g_kineticstress" in name:
                val = self.global_obs['g_kineticstress'].getVal()[int(name[-2]),int(name[-1])]
            elif "g_borntensor" in name:
                val = self.global_obs['g_borntensor'].getVal()[int(name[-2]),int(name[-1])]
            self.file.write("{:.12f} ".format(val/self.sysdata.simbox.vol))
        self.file.write("\n")
    
    def update(self,frame_num):
        self.sysdata.update(frame_num);
        #Update the global calculator
        self.ikglobal = _pyglasstools.GlobalCalculator(self.sysdata._getParticleSystem(),self.pair._getPairPotential())
        #Then, add observables
        self.global_obs = observables.initialize_global(self.__glob_obs_names)
        for name in self.global_obs:
            self.ikglobal.addObservable(self.global_obs[name]._getObservable())

def analyze(frame_list):
    for frame_num in frame_list:
        for logger in loggers:
            logger.update(frame_num);
            logger.run();
            logger.save(frame_num);
