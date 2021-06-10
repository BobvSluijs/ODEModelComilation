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
from Scoring import *
import statsmodels.api as sm
from PlotData import *

from matplotlib.tri import Triangulation
from scipy.spatial import ConvexHull
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt 
from matplotlib import colors as mcolors
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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
    def __init__(self,model,include = [],conditions = {}, samples = 100,time = (0,800),glb = '',dt = 1,sf = "oscillations",name = '',likelihood = False,control_space = False,plot =True):
        """set which parameters are taken into account in this analysis"""
        if not include:
            include = model.include
        if control_space:
            include = model.control_parameters
        include = model.control
#        include = list(model.p_space.keys())
        if 'kflow' in include:
            include.remove('kflow')


#!!	N_CymR	N_PhlF	N_TetR	S70	degmRNACymR	degmRNAPhlF	degmRNAS19	degmRNATetR	degmRNAdeGFP	degmRNAmmCherry	dil	kd_CymR	kd_PhlF	kd_S19	kd_S70	kd_TetR	kd_p19CymR	kd_p19PhlF	kd_p19TetR	kmatdeGFPdark	kmatmmCherry
            
        V ={'KcatP19':1 ,
        'KcatP70':1.2 ,	
        'KcatmRNACymR': 12,
        'KcatmRNAPhlF':9,	
        'KcatmRNAS19': 4,
        'KcatmRNATetR': 12,	
        'KcatmRNAdeGFP': 4,			
        'N_CymR':4.5 ,	
        'N_PhlF': 1.2,	
        'N_TetR': 2.3,	
        'S70':43 ,
        'degmRNACymR':0.1 ,	
        'degmRNAPhlF': 0.1,	
        'degmRNAS19': 0.014,	
        'degmRNATetR': 0.2,
        'degmRNAdeGFP': 0.07,
        'degmRNAmmCherry': 0.6,	
        'dil': 0.02,	
        'kd_CymR':40,	
        'kd_PhlF': 200,	
        'kd_S19':4 ,	
        'kd_S70':50 ,	
        'kd_TetR': 100,	
        'kmatdeGFPdark':0.09}
        
        for k,v in model.fixed.items():
            if k in V.keys():
                model.fixed[k]= V[k]
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
        self.data    = {} 
        '''The bit of code that iterates through the generated paramater space''' 
        i = 0
        for var in self.glb:
            print(i)
            """create the parameter vectors and initial corditions to solve the system"""
            variables = ModelVariables(model,modification = var,conditions = conditions)    
            """"solve the system"""
            solution = ModelSolver(model,variables = variables,simtime = time,dt = dt,forward_sensitivity = False)
            """add the sensitivities to the vector"""
            data = solution.__getData__()    
            """run it again to get bistability with previous conditions"""
            self.data[i] = assessment['oscillations'](data.simdata,observables = model.observables)
            if plot:
#                count = len(os.listdir('C:\\Users\\Bob\\Desktop\\Oscillations 2\\'))
                try:
                    os.makedirs('C:\\Users\\Bob\\Desktop\\Demo Simulation\\')
                    print('success')
                except:
                    pass
                count = len(os.listdir('C:\\Users\\Bob\\Desktop\\Demo Simulation\\'))
                for ID,timeseries in data.simdata.items():
                    plt.plot(timeseries,label = ID)
                    plt.legend(fancybox = True)
                    plt.xlabel("Time")
                    plt.ylabel("Concentration")
                plt.savefig('C:\\Users\\Bob\\Desktop\\Demo Simulation\\{}.png'.format(str(count)))
                plt.close()
            """update counter and see it fly"""
            i += 1
        


    def analyse(self,model,store = True,distributions = True, statistics = True, covariance = True, phase = False,plot =[]):
        '''every scoring function has a double return function which returns the criteria of oscillations including
        scores which have a binary proposition i.e. oscillations yes or no''' 
        criteria,binary = assessment[self.sf]([],[],criteria = True)
        '''function in transform data to get the a somewhat more wieldy format to obtain easily plotted arrays''' 
        p_all,s_all,p_bin,s_bin = unpack(self.glb,self.data,model.observables,criteria = criteria,binary = binary)

        if not plot:
            plot = p_all.keys()
        """ get the quantiles for the datasets to assess correlations in a PCA analysis """ 
        if distributions:
            ''' Loop through the observables and the criteria '''
            for st in model.observables:
                if len(criteria) > 2 and len(criteria)%2 == 0:
                    clm = len(criteria)/2
                else:
                    clm = len(criteria)
    
                fig = plt.figure(figsize=(10,10))  
                row = 1; 
                ''' similar to the local sensitivity we create lists of the criteria and fill up the scores '''
                info = {}
                for i in criteria:
                    info[i] = []
                for i in range(len(self.data)):
                    for k,v in self.data[i][st].items():
                        info[k].append(v)   
                        
                """ The plots are different for binary criteria and for quantitative measurements
                we start with the quanititive criteria"""
                for i in criteria:
                    if i not in binary:
                        ax = plt.subplot(len(criteria),clm,row)
                        bl = max(info[i]) - min(info[i]) + 10
                        if bl >= 100:
                            bl = 50
                        ax.hist(info[i],bins = int(bl),color = 'tab:blue',edgecolor = 'k',alpha = 0.8)
                        ax.set_title("{0} \n State: {1}".format(i,st),size = 14)
                        row += 1
                        
                """ next we plot the binary criteria"""                    
                for i in criteria:
                    if i in binary:
                        ax = plt.subplot(len(criteria),clm,row)
                        pos = [j for j in info[i] if j == True]
                        neg = [j for j in info[i] if j == False]
                        total    = len(pos) + len(neg)
                        ax.bar([0.5,1.5], [100*float(len(pos))/float(total),100*float(len(neg))/float(total)], 
                               align='center',color = 'tab:blue',edgecolor = 'k',alpha = 0.8)
                        ax.set_xlabel("   +                               -   ",size = 14)   
                        ax.set_ylabel("Percentage %",size = 14)
                        ax.set_title("{0} \n State: {1}".format(i,st),size = 14)
                        row += 1
                row = 0
                
                '''Save the figure to the file'''
                if store:
                    plt.tight_layout()
                    plt.savefig(model.figure_path + "global_distribution_" +  filecount(model.figure_folder),dpi = 600)
#                plt.show()   
                
        if statistics:    
            '''Save the figure to the file'''        
            critdict = {i:[] for i in criteria}; bindict  = {i:[] for i in binary}
            """we create a set of dictionaries to get the right correlations between parameters 
            specifically we loop through the data and we take the criteria that correspond to 
            the True in the binary data and quantitative data i.e. scoring function that have a gradient"""
            correlation = {i:copy.deepcopy(critdict) for i in model.observables} 
            if binary:                              
                logit_correlation = {i:copy.deepcopy(bindict) for i in model.observables} 
                for state in model.observables:
                #calulate the statistical variation of the items
                    """ we calculate the pearson correlation between scores and parameters """
                    for parameter,var in p_bin[state].items():
                        for criterion,out in s_bin[state].items():
                                var = numpy.array(var)
                                out = numpy.array(out)
                                pearson,p_value = scipy.stats.pearsonr(var, out)
                                correlation[state][criterion].append((parameter,pearson,p_value))
                    """ A binomial fit of the True/False data with the parameters, Note that not all the conditions/assumptions
                    for this model are met and should serve as a rough comparison"""             
                    for parameter,var in p_all.items():
                        for criterion in binary:
                            var = numpy.array(var)
                            out = numpy.array(s_all[state][criterion])
                            glm_binom = sm.GLM(out, var, family=sm.families.Binomial()) 
                            result = glm_binom.fit()
                            logit_correlation[state][criterion].append((parameter,numpy.exp(result.params)[0]))   
            """we write in the exception if there is no binary value for which the data needs to be seperated"""    
                #calulate the statistical variation of the items                    
            if not binary:
                for parameter,var in p_all.items():
                    for state,scores in s_all.items():
                        for criterion,out in scores.items():
                            var = numpy.array(var)
                            out = numpy.array(out)
                            pearson,p_value = scipy.stats.pearsonr(var, out)
                            correlation[state][criterion].append((parameter,pearson,p_value))
        
            ''' Loop through the observables and the criteria and plot the things'''
            for state in model.observables:
                if len(criteria) > 2 and len(criteria)%2 == 0:
                    clm = len(criteria)/2
                else:
                    clm = len(criteria)
                ''' We plot the figures in a figure panel specifically barplots '''
                row = 1
                fig = plt.figure(figsize=(10,10))  
                for criterion in criteria:
                    ax = plt.subplot(len(criteria),clm,row)
                    if criterion not in binary:      
                        pearson = [];pvalues = []; parameters = [] 
                        for i in correlation[state][criterion]:
                            pr,corr,sig = i
                            parameters.append(pr)
                            pearson.append(corr)
                            pvalues.append(sig) 
                        title  = 'Pearson Correlation: {}'.format(criterion)
                        ylabel = 'Pearson \n state - {}'.format(state)
                        ax = barplot_1D(pearson,parameters,title = title,yl = ylabel,ax = ax)
                        row += 1
                    ''' loop trough the binary criteria and create logit odds barplots '''               
                    if criterion in binary:
                        pearson = []; parameters = [] 
                        for i in logit_correlation[state][criterion]:
                            pr,corr = i
                            parameters.append(pr)
                            pearson.append(corr)
                        title = 'Logit Odds: {}'.format(criterion)
                        ylabel= 'Odds Ratio \n state - {}'.format(state)
                        ax = barplot_1D(pearson,parameters,title = title,yl = ylabel,ax = ax)
                        row += 1 
                ''' save the plot'''
                plt.tight_layout()         
                if store:  
                    plt.savefig(model.figure_path + "global_distribution_" +  filecount(model.figure_folder),dpi = 600)
#                plt.show()             

        ''' the phase is assess the shift in quantitative variables between 
        combinations of 2 rates and or look at the size of the binary criteria
        '''     
        if phase:  
            try:
                tri = list(itertools.combinations(self.glb[-1].keys(),3))
                for x_id,y_id,z_id in tri:
                    if x_id in plot and y_id in plot and z_id in plot:
                        fig = plt.figure(figsize=(5*len(model.observables),8))
                        ''' loop through the observables and do a convexhull to obtain outer points and
                        the triangulations to connect them '''                    
                        row = 1
                        for st in model.observables:
                            ax = fig.add_subplot(len(model.observables), len(model.observables), row, projection='3d')
                            X = p_bin[st][x_id];Y = p_bin[st][y_id];Z = p_bin[st][z_id];
                            """axis needs to be rescaled if the scale of the parameters is logarithmic"""
                            xdist = check_distribution(model.boundaries[x_id]);ydist = check_distribution(model.boundaries[y_id]);zdist = check_distribution(model.boundaries[z_id]);
                            ax.set_xlabel(x_id,size = 12);ax.set_ylabel(y_id,size = 12); ax.set_zlabel(z_id,size = 12);
                            if xdist == "log":
                                ax.set_xlabel("$log_{10}$("+x_id+")",fontsize = 12)
                                X = numpy.log10(X)
                            if ydist == "log":
                                ax.set_ylabel("$log_{10}$("+y_id+")",fontsize = 12)
                                Y = numpy.log10(Y)
                            if zdist == "log":
                                ax.set_zlabel("$log_{10}$("+z_id+")",fontsize = 12)
                                Z = numpy.log10(Z)
                            cvx = ConvexHull(numpy.array([X,Y,Z]).T)
                            tri = Triangulation(X,Y, triangles=cvx.simplices)
                            ax.plot_trisurf(tri,Z,cmap = "gist_earth")
                            centroid = [numpy.mean(cvx.points[cvx.vertices,i]) for i in range(3)]
                            ax.tick_params(axis='both', which='minor', labelsize=12)
                            ax.set_title("Phase {0} \n {1}: Volume = {2}".format(st,binary[0],str(round(cvx.volume,2))))
                            plt.draw()  
                            row += 1 
                        ''' save the figure in the model file'''         
                        plt.tight_layout()                 
                        if store:
                            plt.savefig(model.figure_path + "convexhull_" +  filecount(model.figure_folder),dpi = 600)
#                        plt.show()   
            except ValueError:
                print("There are no binary criteria so no convex hull calculations please continue")

            ''' create a list which combines all the parameters and plot the scores for each 
            states and each biset in large plts'''                     
            bisets  = list(itertools.combinations(self.glb[-1].keys(),2))
            for bi in bisets:
                if bi[0] in plot and bi[1] in plot:
                   
                    clm    = len(criteria) - len(binary)
                    fig = plt.figure(figsize=(20,20))
                    row = 1; 
                    for state in model.observables:                
                        for cr in criteria:
                            ax = plt.subplot(len(criteria),clm,row)
                            x,y = bi
                            X = numpy.array(p_all[x]);
                            Y = numpy.array(p_all[y]);
                            Z = numpy.array(s_all[state][cr])
                            ax.set_title("Phase {0} \n {1}".format(state,cr))
                            ax.tick_params(axis='both', which='minor', labelsize=12)
                            ax.set_xlabel(x,fontsize = 14)
                            ax.set_ylabel(y,fontsize = 14)
                            
                            """axis needs to be rescaled if the scale of the parameters is logarithmic"""
                            xdist = check_distribution(model.boundaries[x])
                            ydist = check_distribution(model.boundaries[y])                     
                            if xdist == "log":
                                ax.set_xlabel("$log_{10}$("+x+")",fontsize = 12)
                                X = numpy.log10(X)
                            if ydist == "log":
                                ax.set_ylabel("$log_{10}$("+y+")",fontsize = 12)
                                Y = numpy.log10(Y)
                                
                            cp = ax.tripcolor(X, Y, Z,20,cmap = 'gist_earth',antialiased = True)
                            if cr in binary:
                                cp = ax.tripcolor(X, Y, Z,lw=2,cmap = 'gist_earth',antialiased = True,vmin = 0, vmax = 1)
                            fig.colorbar(cp,ax = ax)
                            row += 1 
                    ''' save the figure '''                                       
                    if store:
                        plt.savefig(model.figure_path + "phaseplot_" +  filecount(model.figure_folder),dpi = 600)
#                    plt.show()    
            
        #rewrite the PCA analysis
        if covariance:
            '''create parameter bisets to get correlation'''
            bisets = {i:0 for i in set(list(itertools.combinations(self.glb[-1].keys(),2)))}     
            '''Correlation is observable specific '''             
            correlation = {i:copy.deepcopy(bisets) for i in model.observables}  
            if binary:
                Biplot(p_bin,model.observables,store = True, folder = model.figure_folder, path = model.figure_path)
                for o in model.observables:
                    '''loop through all combinations of the data (but first normalize the parameters in p_bin)''' 
                    normalized  = {p:numpy.array(data) - numpy.mean(data) for p,data in p_bin[o].items()}
                    for p1,p2 in bisets.keys():
                        c1 = normalized[p1]
                        c2 = normalized[p2]
                        pearson,p_value = scipy.stats.pearsonr(c1, c2) 
                        correlation[o][(p1,p2)] = pearson
                    '''gather the data and plot in a heatmap with correct
                    labels and relation'''
                    for cr in binary:
                        heatmap(correlation[o],normalized.keys(),criterion = cr)
                    quant = quantiles(s_bin,p_all,model.observables,criteria,binary,fraction = 0.1)
            else:
                 quant =quantiles(s_all,p_all,model.observables,criteria,binary,fraction = 0.1)  
            
            ''' Looping over the observables and the criteria  '''
            for qnt in quant:
                normalized = p_all.keys()
                corr = {};row = 1;
                clm = (len(criteria)-len(binary))*len(model.observables)/2
                figure = plt.figure(figsize=(20,20))  
                for o in model.observables:
                    corr[o] = {}
                    for cr in criteria:
                        if cr not in binary:
                            corr[o] = {i:'' for i in bisets.keys()}
                            for p1,p2 in bisets.keys():
                                c1 = qnt[o][cr][p1]
                                c2 = qnt[o][cr][p2]
                                '''calculate pearson coefficient anbd add to the dictionary'''
                                pc,p_value = scipy.stats.pearsonr(c1, c2) 
                                corr[o][(p1,p2)] = pc
                            ax = plt.subplot(len(criteria)-len(binary),clm,row)
                            ax,im = heatpanel(corr[o],normalized,criterion = cr,ax = ax)
                            fig.colorbar(im,ax = ax)
                            row += 1
                plt.tight_layout()
                if store:
                    plt.savefig(model.figure_path + "heatpanel_" +  filecount(model.figure_folder),dpi = 600)
#                plt.show()    

            
    
import Repress as demo
models,control = demo.main() 
"""create the SBML models"""
for i in range(len(models)):
    model = models[i]
    """compile the model first to sbml then to C++"""
    model.SBMLconversion()
    model.PytoCompile()

"""two models in model files"""
model_1 = models[0]

"""do a global sensitivty analysis"""
a = GlobalSensitivity(model,plot = True)
#a = a.analyse(model)   
#    
#   
#
#	1.517117117	1.844144144	7.393693694	9.567567568	2.153153153	4.81981982	42.85285285	0.078318318	0.277417417	0.201301301	0.02015015	0.3	0.3	0.0333003	1.006938631	224.5697996	10	136.5007807	2.340827276	242.3172794	1.524695727	17.50827032	0.099023023	0.050833834
#
#{'KcatP19':2.5 ,
#'KcatP70':3.7 ,	
#'KcatmRNACymR': 9,
#'KcatmRNAPhlF':3 ,	
#'KcatmRNAS19': 14,
#'KcatmRNATetR': 15,	
#'KcatmRNAdeGFP': 5,	
#'KcatmRNAmmCherry': ,	
#'Kcatp19CymR': 1.5,
#'Kcatp19PhlF':1.8 ,	
#'Kcatp19TetR': 7,	
#'N_CymR': 9,	
#'N_PhlF': 2,	
#'N_TetR': 4,	
#'S70':43 ,
#'degmRNACymR':0.07 ,	
#'degmRNAPhlF': 0.27,	
#'degmRNAS19': 0.2,	
#'degmRNATetR': 0.02,
#'degmRNAdeGFP': 0.3,
#'degmRNAmmCherry': 0.3,	
#'dil': 0.02,	
#'kd_CymR':1 ,	
#'kd_PhlF': 224,	
#'kd_S19':10 ,	
#'kd_S70':136 ,	
#'kd_TetR': 3,
#'kd_p19CymR': 242,	
#'kd_p19PhlF': 1.5,	
#'kd_p19TetR': 17,	
#'kmatdeGFPdark':0.09}












