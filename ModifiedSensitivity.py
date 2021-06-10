# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 12:04:41 2021

@author: Bob van Sluijs"""


import os 
import sys 

from SolveSystem import ModelSolver
from Model import ModelVariables
from Operations import *
from Parse import *

import time as timing
import matplotlib.pylab as plt
import libsbml
import importlib
import amici
import os
import math
import copy
import sys
import numpy as numpy
    
"""small datacollection of the dataset"""  
class DataCollection:
    def __init__(self,data,parameters):
        self.data       = data
        self.parameters = parameters
        
    def analyse(self,):
        pass
    
class GlobalSensitivity:
    def __init__(self,model,include = [],conditions = {}, samples = 1000,time = (0,4000),glb = '',dt = 1,sf = "ESS",name = '',likelihood = False,control_space = False,plot =True):
        """set which parameters are taken into account in this analysis"""
        if not include:
            include = model.include
        if control_space:
            include = model.control_parameters
        include = model.control
        if 'kflow' in include:
            include.remove('kflow')
        """we add the model to the analysis, create a data storage dict which can be taken
        up by the plotting functions, we create a global space for each parameter""" 
        if glb == '':
            self.glb = build_global_search(model.p_space,include = include,samples = samples, order = int(math.log10(model.spacing)))
        else:
            self.glb = glb
            
        self.sf = sf
        self.newvars = []
        '''include the parameters to be included in the generation of the space and
        Build a dictionary with the relevant scoring data organized by:
        parameter set in phase, observable and criteria, the latter being in evaluate
        i.e. {i:{state1:{crit1:score,crit2:score etc} etc}'''
        self.data    = []   
        '''The bit of code that iterates through the generated paramater space''' 
        i = 0
        for var in self.glb:
            print(i)
            """create the parameter vectors and initial corditions to solve the system"""
            variables = ModelVariables(model,modification = var,conditions = conditions)    
            """"solve the system"""
            solution = ModelSolver(model,variables = variables,simtime = time,dt = dt)
            """add the sensitivities to the vector"""
            data = solution.__getData__()    
            """run it again to get bistability with previous conditions"""
            self.data.append(DataCollection(data.simdata,var))
            if plot:
                for ID,timeseries in data.simdata.items():
                    plt.plot(timeseries,label = ID)
                    plt.legend(fancybox = True)
                    plt.xlabel("Time")
                    plt.ylabel("Concentration")
                    plt.show()
    

import SubstrateCompetitionModels as MMC
models,control = MMC.main() 
"""create the SBML models"""
for i in range(len(models)):
    model = models[i]
    """compile the model"""
    model.SBMLconversion()
    model.PytoCompile()

GlobalSensitivity(models[0])
    
    
    
    

    













