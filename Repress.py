# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 12:45:09 2021

@author: Bob
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 22:39:33 2021

@author: Bob
"""
from Model import ModelObject
class optimization_information:
    def __init__(self,theta = [],score = [],fitfunc = [],include = [],conditions = [],plot = []):
        self.theta      = theta
        self.score      = score
        self.fitfunc    = fitfunc
        self.include    = include
        self.conditions = conditions
        self.plot       = plot
        
        
        
def main():
    models,control = {},{}
    master_boundaries = {'(in)p70_S19':(2,6),
                         '(in)p19CymR_TetR':(0.1,2),
                         '(in)p19Phlf_CymR':(0.1,2),
                         '(in)p19TetR_PhlF':(1,5),
                         '(in)p19TetR_deGFP':(4,8),
            
            
                        'S70': (40,45),
                        'KcatP70':(0.1,10),
                        'kd_S70':(1,1000),
                        'degmRNAS19':(0.01,0.3),
                        'KcatP19':(0.1,10),
                        'kd_S19':(1,1000),
                        'degmRNAdeGFP':(0.01,0.3),    
                        'degmRNAmmCherry':(0.01,0.3),
                        'KcatmRNAS19':(0.1,10), 
                        'KcatmRNAmmCherry':(0.1,10),
                        'KcatmRNAdeGFP':(0.1,10),
                        'kmatdeGFPdark':(0.099,0.1),
                        'kmatmmCherry':(0.05,0.051),
                       
                        'degmRNACymR':(0.01,0.4),    
                        'kd_CymR':(1,1000),
                        'N_CymR':(1,10),
                        'Kcatp19CymR':(0.1,10),
                        'kd_p19CymR':(1,1000),
                        'KcatmRNACymR':(0.1,10), 
     
                        'degmRNAPhlF':(0.01,0.4),    
                        'kd_PhlF':(1,1000),
                        'N_PhlF':(1,10),
                        'KcatmRNAPhlF':(0.1,10),
                        'Kcatp19PhlF':(0.1,10),
                        'kd_p19PhlF':(1,1000), 
 
                        'degmRNATetR':(0.01,0.4),    
                        'kd_TetR':(1,1000),
                        'N_TetR':(1,10),
                        'KcatmRNAS19':(0.1,10),
                        'KcatmRNATetR':(0.1,10),
                        'Kcatp19TetR':(0.1,10),
                        'kd_p19TetR':(1,1000),
                        'dil':(0.026,0.043)}   
     
    master_fixed   = {'(in)p70_S19':1,
                         '(in)p19CymR_TetR':1,
                         '(in)p19Phlf_CymR':1,
                         '(in)p19TetR_PhlF':1,
                         '(in)p19TetR_deGFP':1,          
                       
                       'S70':40,
                       'KcatP70':1,
                       'kd_S70':50,
                       'degmRNAS19':0.2,
                       'KcatP19':1,
                       'kd_S19':50,
                       'degmRNAdeGFP':0.2,    
                       'degmRNAmmCherry':0.2,
                       'KcatmRNAS19':1, 
                       'KcatmRNAmmCherry':1,
                       'KcatmRNAdeGFP':1,
                       'kmatdeGFPdark': 0.1,
                       'kmatmmCherry':0.05,
                       
                       'degmRNACymR':0.2,    
                       'kd_CymR':50,
                       'N_CymR':2,
                       'Kcatp19CymR':1,
                       'kd_p19CymR':50,
                       'KcatmRNACymR':1, 
    
                       'degmRNAPhlF':0.2,    
                       'kd_PhlF':50,
                       'N_PhlF':2,
                       'KcatmRNAPhlF':1,
                       'Kcatp19PhlF':1,
                       'kd_p19PhlF':50, 

                       'degmRNATetR':0.2,    
                       'kd_TetR':50,
                       'N_TetR':2,
                       'KcatmRNAS19':1,
                       'KcatmRNATetR':1,
                       'Kcatp19TetR':1,
                       'kd_p19TetR':50,
                       'dil':0.026}   
    
    
    parameters = ['kd_PhlF','kd_TetR','kd_CymR','S70','(in)p70_S19','(in)p19CymR_TetR','(in)p19Phlf_CymR','(in)p19TetR_PhlF','(in)p19TetR_deGFP','dil',
                  'KcatP70','kd_S70','kd_S70','degmRNAS19','KcatP19','kd_p19PhlF','kd_p19TetR','kd_p19CymR','N_PhlF','N_TetR','N_CymR','degmRNACymR','degmRNAPhlF',
                  'degmRNATetR','degmRNAdeGFP','KcatmRNAS19','KcatmRNACymR','KcatmRNATetR','KcatmRNAPhlF','KcatmRNAdeGFP','kmatdeGFPdark']
    fixed      = {i:master_fixed[i] for i in parameters}
    boundaries = {i:master_boundaries[i] for i in parameters}
    modelname = 'ReconstructedRepressilatorTetReGFPfour'  
    stringmodel = """ + (in)p70_S19 * dil - p70_S19 * dil 
                      + (in)p19CymR_TetR * dil - p19CymR_TetR * dil 
                      + (in)p19Phlf_CymR * dil - p19Phlf_CymR * dil 
                      + (in)p19TetR_PhlF * dil - p19TetR_PhlF * dil  
                      + (in)p19TetR_deGFP * dil - p19TetR_deGFP * dil 
                      + p70_S19 * KcatP70 * ( S70 / kd_S70 ) / ( 1 |+| ( S70 / kd_S70 ) ) - dil * mRNAS19 - degmRNAS19 * mRNAS19 
                      + p19Phlf_CymR * KcatP19 * ( S19 / kd_p19PhlF ) / ( 1 |+| ( S19 / kd_p19PhlF ) |+| ( PhlF / kd_PhlF ) ** N_PhlF )  - dil * mRNACymR - degmRNACymR * mRNACymR  
                      + p19TetR_PhlF * KcatP19 * ( S19 / kd_p19TetR ) / ( 1 |+| ( S19 / kd_p19TetR ) |+| ( TetR / kd_TetR ) ** N_TetR )  - dil * mRNAPhlF - degmRNAPhlF * mRNAPhlF 
                      + p19CymR_TetR * KcatP19 * ( S19 / kd_p19CymR ) / ( 1 |+| ( S19 / kd_p19CymR ) |+| ( CymR / kd_CymR ) ** N_CymR )  - dil * mRNATetR - degmRNATetR * mRNATetR          
                      + p19TetR_deGFP * KcatP19 * ( S19 / kd_p19TetR ) / ( 1 |+| ( S19 / kd_p19TetR ) |+| ( TetR / kd_TetR ) ** N_TetR ) - dil * mRNAdeGFP - degmRNAdeGFP * mRNAdeGFP                       
                      + KcatmRNAS19 * mRNAS19 - dil * S19 
                      + KcatmRNACymR * mRNACymR - dil * CymR 
                      + KcatmRNAPhlF * mRNAPhlF - dil * PhlF 
                      + KcatmRNATetR * mRNATetR - dil * TetR    
                      + KcatmRNAdeGFP * mRNAdeGFP - kmatdeGFPdark * deGFPdark - dil * deGFPdark 
                      + kmatdeGFPdark * deGFPdark - dil * deGFP    """
                      
    states             = ['p70_S19','p19Phlf_CymR','p19CymR_TetR','p19TetR_PhlF','p19TetR_deGFP','mRNAS19',
                          'mRNACymR','mRNAPhlF','mRNATetR','mRNAdeGFP',
                          'S19','CymR','PhlF','TetR','deGFPdark','deGFP']
    observables        = ['deGFP']
    control_parameters = ['(in)p70_S19','(in)p19CymR_TetR','(in)p19Phlf_CymR','(in)p19TetR_PhlF','(in)p19TetR_deGFP','dil']
    conditions         = [{'(in)p70_S19':1,'(in)p19CymR_TetR':1,'(in)p19Phlf_CymR':1,'(in)p19TetR_PhlF':1,'(in)p19TetR_deGFP':1,'dil':0.026}]
    models[len(models)]     = ModelObject(stringmodel,states,boundaries,fixed,observables = observables,name = modelname,control_parameters = control_parameters)   
    control[len(control)]   = optimization_information(conditions = conditions)  
    return models,control
