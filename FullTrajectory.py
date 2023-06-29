import tellurium as te
import roadrunner
import libsbml
from collections import OrderedDict
import numpy as np
import os
import pickle as pk
import csv
from PyQt5.QtWidgets import QFileDialog
#from roadrunner import Config as rrc
#rrc=MAX_OUTPUT_ROWS=1000000000

class Gil_Trajectory():
    def __init_(self):
        self.systemdict=''
        self.savehistoarray=''
        self.model=''
        
        self.allreactions=''
        self.allstates=''
        self.histoarrays=''
        self.modelspecies=[]

    #def model(self):
        #pass
            
    def trajectory(self,duration=2000):
        '''creates a trajectory of events and returns an array of the population of each species and elapsed time using tellurum'''
        te.setDefaultPlottingEngine('matplotlib')
        #sbml_model = te.antimonyToSBML(self.model)
        #r=te.loada(self.model)
        self.savehistoarray = os.path.dirname(os.path.abspath(__file__))+'\\histoarrays'
        if not os.path.exists(self.savehistoarray):
            os.mkdir(self.savehistoarray)
        try:
            r = te.loada(self.model)
        except:
            r = te.loadSBMLModel(self.model)
        r.setIntegrator('gillespie')
        r.integrator.setValue('max_output_rows',100000000)
        r.integrator.setValue('nonnegative', True)
        r.integrator.seed=147258369
        self.modelspecies=r.selections[1:]
        results=[]
        results=r.gillespie(0, duration)
        # for n,t in enumerate(results1):
        #     if t[0]>2:
        #         results=results1[n:]
        #         break
        allreactions=[]
        allstates=[]
        system1=OrderedDict()
        oldTime=results[0][0]
        oldState=[int(x) for x in results[0][1:]]
        for i in results[1:]:
            reaction=[int(x) for x in np.array(i[1:])]-np.array(oldState)
            if not all(v == 0 for v in reaction):
                timestep=i[0]-oldTime
        
                if not reaction.tolist() in allreactions:
                    allreactions.append(reaction.tolist()) ## appends a reaction vector in the list if it isn't already in the list. 
        
                if not oldState in allstates:
                    allstates.append(oldState)
                if not [int(x) for x in i[1:]] in allstates:
                    allstates.append([int(x) for x in i[1:]])

                if tuple(oldState) in system1.keys():  
                    system1[tuple(oldState)].append([timestep, allstates.index([int(x) for x in i[1:]])])  
                else:  
                    system1[tuple(oldState)] = [[timestep,allstates.index([int(x) for x in i[1:]])]]
                    
                oldState =[int(x) for x in i[1:]]
                oldTime = i[0]
            else:
                oldState =[int(x) for x in i[1:]]
                oldTime = i[0]
        self.systemdict=system1
        self.allreactions=allreactions
        self.allstates=allstates
        self.histoarrays=self.histoarray() 
        Results = os.path.dirname(os.path.abspath(__file__))+'\\Results\\'
        if not os.path.exists(Results):
            os.mkdir(Results)
        # fullpath2=os.path.join(Results,'Allstates.txt')
        # with open (fullpath2.csv,'w') as f23:
        #     wr=csv.writer(f23)
        #     wr.writerows(self.allstates)    
    def histoarray(self):     
        ''' creates an array of the elapsed time of the observed events and the number of observed events'''           
        histoarrays=[]
        for k,v in self.systemdict.items():
            histoarray=v
            histoarrays.append([k,histoarray])
            savefile1="histoarray{}.db".format(k)
            fullpath1=os.path.join(self.savehistoarray,savefile1)
            with open(fullpath1, "ab") as datahis:
                pk.dump([k,histoarray], datahis)
        return histoarrays
        
if __name__ == "__main__":
    tr=Gil_Trajectory()
    tr.trajectory()
    

