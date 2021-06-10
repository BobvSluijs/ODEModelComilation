# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 13:39:57 2019

@author: Bob van Sluijs
"""

from Model import ModelObject
"""A Super Main file meant to test a lot of individial models versus the pipeline that 
was developed, a benchmark to spot inconsistancies essentially"""
class optimization_information:
    def __init__(self,theta = [],score = [],fitfunc = [],include = [],conditions = [],plot = []):
        self.theta      = theta
        self.score      = score
        self.fitfunc    = fitfunc
        self.include    = include
        self.conditions = conditions
        self.plot       = plot
    
        
        
def main():   
    models = {}
    control = {}
    ##### you build a function that defines a space for a givenset of bounfaries 
    description = "{} this model analysis was performed by bob to see ascertain the oscillating region in a dual reactor environment"
    #Define parameter boundaries
    boundaries = {"KmTrCg" :(0.0706,7.06), "KmTrS1":(0.0738,7.06), "KcatTrS1" :(128,12800),
    "KmCrS1" : (0.0913,0.0913),"KcatTrCg": (29,29),"KcatCrS1" : (735,735),"kflow" : (0.01,0.011),
    "Tr0" : (0.0000001,0.001),"Cg0" : (0.0005,5), "S10" : (5,100)}
    
    fixed = {"KmTrCg" :0.706,"KmTrS1": 0.738,"KcatTrS1" : 1280,"KmCrS1" : 0.0913,
    "KcatTrCg": 29,"KcatCrS1" : 735,"kflow" : 0.01,"Tr0" : 0.001,"Cg0" : 0.001,"S10" : 5}
    #define control parameters
    
    control_parameters = ['Tr0','Cg0','S10','kflow'] #what are the control parameters
    #define equations of you syste  m
    
    stringmodel = ''' + kflow * Tr0 - kflow * Tr 
    - ( KcatTrCg * Tr * Cg ) / ( KmTrCg * ( 1 |+| ( S1 / KmTrS1 ) ) |+| Cg )  + kflow * Cg0 - kflow * Cg 
    + ( KcatTrCg * Tr * Cg ) / ( KmTrCg * ( 1 |+| ( S1 / KmTrS1 ) ) |+| Cg )  - kflow * Cr 
    - ( KcatTrS1 * Tr * S1 ) / ( KmTrS1 * ( 1 |+| ( Cg / KmTrCg ) ) |+| S1 )  - ( KcatCrS1 * Cr * S1 ) / ( KmCrS1 |+| S1 ) + kflow * S10 - kflow * S1 
    + ( KcatCrS1 * Cr * S1 ) / ( KmCrS1  |+| S1 ) - kflow * Or 
    + ( KcatTrS1 * Tr * S1 ) / ( KmTrS1 * ( 1 |+| ( Cg / KmTrCg ) ) |+| S1 ) - AMC * kflow '''
    
    #define states of equations, make sure they are in the same order (or that you define them seperately in a dictionary with equations)
    states      = ['Tr','Cg','Cr','S1','Or','AMC']
    observables = ['AMC']
    #define the scoring function you wish to use
    conditions            = [{"kflow":0.003,"Tr0":0.005,"Cg0":0.0,"S10":0.5}]
    modelname             = 'BistabilityCompetition'
    include               = []
    models[len(models)]   = ModelObject(stringmodel,states,boundaries,fixed,observables = observables,name = modelname,control_parameters = control_parameters)   
    control[len(control)] = optimization_information(conditions = conditions)      #############################################################################################################################################################    
    #############################################################################################################################################################
    #############################################################################################################################################################
    #############################################################################################################################################################
    ##### you build a function that defines a space for a givenset of bounfaries 
    description = "{} this model analysis was performed by bob to see ascertain the oscillating region in a dual reactor environment"
    
    fixed = {"KmTrCg" :0.706,"KmTrS1": 0.738,"KcatTrS1" : 1280,"KmCrS1" : 0.0913,
    "KcatTrCg": 29,"KcatCrS1" : 735,"kflow" : 0.05,"Tr0" : 0.001,"Cg0" : 0.001,"S10" : 5}
    #define control parameters
    
    control_parameters = ['Tr0','Cg0','S10','kflow'] #what are the control parameters
    #define equations of you syste  m
    
    stringmodel = ''' + kflow * Tr0 - kflow * Tr 
    - (  KcatTrCg * Tr * Cg ) / ( KmTrCg |+| Cg )  + kflow * Cg0 - kflow * Cg 
    + ( KcatTrCg * Tr * Cg ) / ( KmTrCg |+| Cg )  - kflow * Cr 
    - ( KcatTrS1 * Tr * S1 ) / ( KmTrS1 |+| S1 )  -  ( KcatCrS1 * Cr * S1 ) / ( KmCrS1 |+| S1 )  + kflow * S10 - kflow * S1 
    + ( KcatCrS1 * Cr * S1 ) / ( KmCrS1  |+| S1 ) - kflow * Or 
    + ( KcatTrS1 * Tr * S1 ) / ( KmTrS1 |+| S1 )  - AMC * kflow '''
    
    #define states of equations, make sure they are in the same order (or that you define them seperately in a dictionary with equations)
    states              = ['Tr','Cg','Cr','S1','Or','AMC']
    observables         = ['AMC']
    #define the scoring function you wish to use
    conditions          = [{"kflow":0.003,"Tr0":0.005,"Cg0":0.0,"S10":0.5}]
    modelname           = 'BistabilityNoCompetition' 
    include             = []
    models[len(models)] = ModelObject(stringmodel,states,boundaries,fixed,observables = observables,name = modelname,control_parameters = control_parameters)   
    control[len(control)] = optimization_information(conditions = conditions)      
    return models,control
