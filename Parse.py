# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:53:31 2019

@author: huckg
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 06 14:11:45 2019

@author: huckg
"""



import numpy 
import os
import csv

from collections import defaultdict
from DataTransform import *
from pyDOE import *


#*__________________excell format for prebiotic reaction networks______________*
class ExcellNetwork:
    def __init__(self,sequence = []):
        """sequence of reactions used to generate smarts"""
        self.sequence   = sequence
        """reactions/secies smiles in the network"""
        self.reactions  = [] 
        """"reaction type e.g. enolate formation, hydration etc"""
        self.rtype      = []
        """if the rate types are the same they will be classified as a single rate governing all reactions of that type"""
        self.itype      = []
        """"forward and or reverse rates for the reaction networks"""
        self.rates = []
        """boundaries"""
        self.boundaries = []
        
    def update(self,rct,rtype,bnd,rti):
        self.reactions.append(rct)
        self.rtype.append(rtype)
        self.itype.append(rti)    
        self.rates.append(bnd)

        
def read_csv_smarts(fname):
    """list of smart definitions"""
    smarts = []
    """open folder and read the file"""
    with open(fname,'r') as csvfile:
        reader = csv.DictReader(csvfile,delimiter=';')
        for row in reader:
            m = (row['SMARTS'],row['Conventional'],row['Database'])
            smarts.append(m)
    compatibility = {i:j for n,i,j in smarts}
    smarts = {i:n for n,i,j in smarts}
    return compatibility,smarts
    
def find_paths(folder):
    return [folder + "\\" + name for name in os.listdir(folder) if os.path.isdir(folder)]

""""this script works with prebiotic reaction networks and asesses the predictive power of said networks
we will, however use excell files with the reaction equations to deal with these types of reactions"""
#*________________Parse the Networks_________________*    
def read_csv_network(folder):
    from NetworkSpace import ReactionNetwork
    """the paths"""
    paths =  [folder + i for i in os.listdir(folder)]
    """the name of the system"""
    names =  [i for i in os.listdir(folder)]
    # Read the reactions from the csv file
    """store networks in this dictionary"""
    networks = {}
    
    """loop through the file with all the networks"""
    count = 0 #topology count
    for path in paths:
        """this tests where the sequence is at the top
        of the input excell sheet with reaction equations"""
        sequence = None
        with open(path,'r') as csvfile:
            data = [line.split() for line in csvfile] 
            for i,x in enumerate(data):
                if i == 0:
                    sequence = x

        network = []
        if "SR1" in sequence[0]:
            """parse the networks generated by the network generator"""            
            with open(path,'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    
                    """set of bools to deconstruct the compounds"""
                    forward_reaction = ''
                    SR1 = row["SR1"]
                    if SR1:
                        forward_reaction += SR1
                    SR2 = row["SR2"]
                    if SR2:
                        forward_reaction += '.' + SR2
                    SR3 = row["SR3"]
                    if SR3:
                        forward_reaction += '.' + SR3
                    forward_reaction += '>>'
                    SP1 = row["SP1"]
                    if SP1:
                        forward_reaction +=  SP1
                    SP2 = row["SP2"]
                    if SP2:
                        forward_reaction += '.' + SP2
                    SP3 = row["SP3"]
                    if SP3:
                        forward_reaction += '.' + SP3
                        
                    """set of reaction types"""                    
                    reaction_type = row["Reaction Type"]
                    """reaction rate"""
                    forward_rate = row["Kf"]
                    """reverse reaction"""
                    reverse_rate = row["Kr"]
                    """reverse reaction"""
                    splt = forward_reaction.split(">>")
                    """reverse reaction"""
                    reverse_reaction = splt[1] + ">>" + splt[0]
                    
                    """network update with the forward and revesrse reaction
                    forward:"""
                    network.append((forward_reaction,reaction_type,forward_rate,''))
                    """reverse:"""
                    network.append((reverse_reaction,reaction_type,reverse_rate,''))

        else:
            """parse the networks generated by the network generator"""            
            with open(path,'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    """the forward reactions"""
                    forward_reaction,rtype,lower,upper = (row['reaction'],row['type'], row['lower rate limit/ s-1 or M-1s-1'],row['upper rate limit/ s-1 or M-1s-1'])
                    """update the network, only the forward in this instance"""
                    network.append((forward_reaction,rtype,'',(lower,upper)))
                             
        """update topology count"""
        networks[count] = network   
        count += 1
        
    """The names of the networks"""
    names    = [i.split(".")[0] for i in names]
    """The network objects with th information"""
    networks = {names[k]:ReactionNetwork(v,name  = names[k]) for k,v in networks.items()}    
    print(networks)
    return networks
    
#*________________Parse the Networks_________________*    
def read_csv_network_generated(folder):
    fnames = find_paths(folder)
    # Read the reactions from the csv file
    """store networks in this dictionary"""
    networks = {}
    
    """loop through the file with all the networks"""
    count = 0 #topology count
    for fname in fnames:
        sequence = None
        with open(fname,'r') as csvfile:
            data = [line.split() for line in csvfile] 
            for i,x in enumerate(data):
                if i == 1:
                    sequence = x
                    
        """parse the networks generated by the network generator"""            
        network = __excellnetwork__(sequence = sequence)
        with open(fname,'r') as csvfile:
            for i in range(2):
                csvfile.next()
            reader = csv.DictReader(csvfile)
            for row in reader:
                rct,rtype,lower,upper,rti = (row['reaction'],row['type'], row['lower rate limit/ s-1 or M-1s-1'],row['upper rate limit/ s-1 or M-1s-1'],row["assignment"])
                network.update(rct,rtype,(lower,upper),rti)
        networks[count] = network    
        """update topology count"""
        count += 1
    return networks

def read_csv_file(fname):
    import os
    # Read the reactions from the csv file
    path = os.path.expanduser(fname)
    reactions_init = []
    with open(fname,'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
          
             reactants = [row['SR1'],row['SR2'], row['SR3']]
             products = [row['SP1'], row['SP2'], row['SP3']]
             reactants = [ r for r in reactants if r != '']
             products = [ p for p in products if p != '']
             if reactants and products:
                 reactions_init.append([reactants, products, row['Kf']])
                 reactions_init.append([products, reactants, row['Kr']])

    # Define species dictionary
    species = {}
    species_tmp = []
    for r in reactions_init:
        for react in r[0]:
            if react in species_tmp:
                continue
            else:
                species_tmp.append(react)
        for prod in r[1]:
            if prod in species_tmp:
                continue
            else:
                species_tmp.append(prod)

    for i in range(len(species_tmp)):
        species[species_tmp[i]] = i

    return reactions_init, species

#*________________Parse the ExperimentalData_________________*    
def ParseExperimentalData_cvs(path,limit = 1e100,scale = '',fitfunc = "leastsquaresfit",delimiter = ','):
    from Measurements import MeasurementObject
    '''Gets a set of datasets from files. Creates dictionary, structure:
    {'Dataset': { 'Conditions':[ [] ], 'Data':{'variable': np.array([])} }  }'''
    data_container = {}
    n = 0
    with open(path,'r') as f:
        for line in f:
            if 'Dataset' in line:
                d_set = line.split(delimiter)[1]
                data_container[d_set] = {'Conditions':[], 'Data':[]}
            n += 1
            if n == limit:
                break
            
    """loop through it the first time and assess conditions"""
    n = 0
    with open(path,'r') as f:
        condition = False
        sets = list(data_container.keys())
        s = 0
        for line in f:
            if 'start_conditions' in line:
                condition = True
                line = next(f)
            if 'end_conditions' in line:
                data_container[sets[s]]['Conditions'] = {k:v for k,v in zip(data_container[sets[s]]['Conditions'][0],data_container[sets[s]]['Conditions'][1])}
                condition = False
                s+=1
            if condition:
                inp = line.strip('\n')
                add = [g for g in inp.split(delimiter) if g != '']
                data_container[sets[s]]['Conditions'].append(add)
            n += 1
            if n == limit:
                break


    """loop through it second time for datavectors """
    n = 0
    with open(path,'r') as f:
        condition = False
        sets = list(data_container.keys())
        s = 0
        for line in f:
            if 'start_data' in line:
                condition = True
                line = next(f)

            if 'end_data' in line:
                condition = False
                s+=1
            if condition == True:
                inp = line.strip('\n')
                add = [g for g in inp.split(delimiter) if g != '']
                data_container[sets[s]]['Data'].append(add)
            n += 1
            if n == limit:
                break
            
    """store both"""
    for d in data_container.keys():
        arr = data_container[d]['Data']
        data = [list(i) for i in zip(*arr[1:])]
        temp_dict = {k:v for k,v in zip(arr[0],data)}
        for k in temp_dict:
            for z in temp_dict[k]:
                new = numpy.array([float(i) for i in temp_dict[k] if i != ''])
                temp_dict[k] = new
        data_container[d]['Data'] = temp_dict
    """parse the datacontainer from the csv file into a measurement class with the 
    appropriate item calls"""
    measurements = {}
    count = 0
    for name,data in data_container.items():
        rawdata = {};units = {}
        for cnd,vector in data["Data"].items():
            if "time" in cnd:
                time = vector
                time_unit = cnd.split(" ")[-1].strip()
            if "M" in cnd:
                smiles = cnd.split(" ")[0]
                rawdata[smiles] = vector
                units[smiles] = [""+cnd.split(" ")[i] for i in range(len(cnd.split(" "))) if i != 0][0]
                
        """make sure we obtain the conditions of the measurement"""
        conditions = {}
        temperature = None; flowrate = None;
        for smiles,cnd in data["Conditions"].items():
            if "flow" in smiles:
                flowrate = eval(cnd)
            elif "Flow" in smiles:
                flowrate = eval(cnd)
            elif "Temperature" in smiles:
                temperature = eval(cnd)
            elif "Temperature" not in smiles and "flow" not in smiles and "time" not in smiles:
                conditions[smiles.split(" ")[0].strip()] = eval(cnd)
 
        """dictionary with measurement object contianing all the experimental data AND settings to be used by the algorithm"""
        m = MeasurementObject(rawdata,time,name = name, conditions = conditions,
                            temperature = temperature,flowrate= flowrate,time_unit = time_unit,units = units,
                            scale = scale, fitfunc = fitfunc,store = True)

        measurements[count] = m
        count += 1      
    """"the measurements are returned"""    
    return measurements


#*________________Optional Experimental Conditions_________________*    
def read_initiators_cvs(path):
    """list of inputs available to intiate the formose reaction"""
    inputs = {}
    """open folder and read the file"""
    with open(path,'r') as csvfile:
        reader = csv.DictReader(csvfile,delimiter=';')
        for row in reader:
            smiles,lower,upper = (row['Species'],row['Lower'],row['Upper'])
            inputs[smiles] = (float(lower),float(upper))
    return inputs


def find_elements(s, ch):
    return [(i,ch) for i, ltr in enumerate(s) if ltr == ch and s[i+1] == ' ']

def define_ratelaw(s):
    ind = list(sorted(find_elements(s,"+")+find_elements(s,"-")))
    rate_equations = []
    for i in range(len(ind)-1):
        fi,fs = ind[i]    #the index of the first sign + or -
        ni,ns = ind[i+1]  #the index of the second sign + or -
        """stringcut which defines the rate"""
        rate_equations.append(s[fi:ni])
    """the final cut"""
    li,ls = ind[-1]
    """"the rate equations to strip"""
    rate_equations.append(s[li:len(s)])
    return rate_equations

#_________________Parse a model to antimony and SBML________________*
def antimony_parser(states,parameters,fluxes,S,ratelaws,staterate,fixed,mapping):
    """find all indices"""
    row,column = numpy.nonzero(S)
    """the paired indices of the S matrix"""
    stoch = [(row[i],column[i]) for i in range(len(row))]
    """find all paired fluxes between states etc."""
    connections = {i:[] for i in range(len(fluxes))}
    for i,j in stoch:
        sign = S[i,j]
        connections[j].append((sign,i)) #state to flux
    fluxlist = []
    """single incoming and outgoing reaction"""
    for flux,inv in connections.items():
        for sign,state in inv:
            reaction = ""
            if sign == -1:
                reaction += " " + states[state] + " -> ; " + fluxes[flux]
            elif sign == 1:
                reaction += " -> {} ; ".format(states[state]) + fluxes[flux]
            fluxlist.append(reaction)
    fluxlist = [i.replace("  "," ") for i in fluxlist]
            
    """create a model list"""
    model = ''
    for i in range(len(fluxlist)):
        """add flux equations"""
        model += "J{}: ".format(i) + fluxlist[i] + ' ; \n'
    model += ''
    """add initial conditions"""
    for i in states:
        model += " "+i+" = 0 ; \n"
    model += ''
    
    """add parameters to the system"""
    for k,v in fixed.items():
        model += " " + k + ' = {} ; \n'.format(str(v))
        
    """the y and p vector replacing the states"""
    y = {states[i]:'y{}'.format(mapping[states[i]]) for i in range(len(states))}
    p = {parameters[i]:'k{}'.format(mapping[parameters[i]]) for i in range(len(parameters))}
    
    """overwrite the states and parameters"""
    for k,v in y.items():
        k = " " + k + " "
        model = model.replace(k,v)
    for k,v in p.items():
        k = " " + k + " "
        model = model.replace(k,v)
        
    """replace the delimiters"""
    model = model.replace("**","^")
    model = model.replace("|+|","+")
    model = model.replace("|-|","-")
    return model