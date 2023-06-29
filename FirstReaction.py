import numpy as np
import tellurium as te
import matplotlib.pyplot as plt
import os
import pickle as pk
import csv
#from roadrunner import Config as rrc
#rrc=MAX_OUTPUT_ROWS=1000000000

class Firstreaction:
    te.setDefaultPlottingEngine('matplotlib')
    def __init__(self, model,mcsteps): 
        self.savehistoarray = os.path.join(os.path.dirname(__file__),'histoarrays')
        if not os.path.exists(self.savehistoarray):
            os.mkdir(self.savehistoarray)
        self.saveResults = os.path.dirname(__file__)+'\\Results'
        if not os.path.exists(self.saveResults):
            os.mkdir(self.saveResults)
        self.systemdict=''
        self.model=model
        self.mcsteps=mcsteps
        #self.endTime=endTime
        self.ResultsDict=''
        self.currentstate=''
        self.allreactions=''
        self.allstates=''
        self.histoarrays=''
        self.probabilities=''

    def ModelResults(self):
        try:
            r = te.loada(self.model)
        except:
            r = te.loadSBMLModel(self.model)
        r.integrator = 'gillespie'
        r.integrator.setValue('max_output_rows',100000000)
        r.integrator.setValue('nonnegative', True)
        r.integrator.seed = 134
        results = []
        for k in range(1, self.mcsteps):
            r.reset()
            s = r.simulate(0, points = 2)#self.endTime,
            results.append(s)
            
            #r.plot(s, alpha=0.7)
        currentState=[int(x) for x in results[0][0][1:]]
        
        return results,currentState

    def first_reaction(self): 
        '''returns a dictionary with keys the states the system visited and entries a list of time elapses for each first reaction''' 
        
        ResultsDict={}
        system1={}
        modelresults=self.ModelResults()
        self.currentstate=tuple(modelresults[1])
        allstates=[list(self.currentstate)]
        system1={}
        for item in modelresults[0]:
            state=tuple(item[1][1:])
            if not (state) in allstates:
                allstates.append(state)
            if state in ResultsDict.keys():
                ResultsDict[state].append(item[1][0]) 
            else:
                ResultsDict[state]=[]
            if self.currentstate in system1.keys():
                system1[self.currentstate].append([item[1][0],state])
            else:
                system1[self.currentstate]=[[item[1][0],state]]
        with open(os.path.dirname(__file__)+'\\allstates.db', "wb") as fostates:
            pk.dump(allstates, fostates)        
        self.ResultsDict = ResultsDict
        self.systemdict=system1
        self.histoarrays=self.histoarray()
        self.allstates=allstates
        allreactions=[]
        for val,i in enumerate(allstates[1:],start=1):
            reactl=np.array(allstates[val])-np.array(allstates[0])
            reactlist=[int(i) for i in reactl]
            allreactions.append(reactlist)

        self.allreactions=allreactions
        self.CondProb()
    
    def CondProb(self):
        '''returns a dictionary with keys the states the system it arrives to and values the probability to go to that state''' 
        ResultsDict=self.ResultsDict
        Cond={}
        Denomin = sum(len(v) for v in ResultsDict.values())
        for state,Nomin in ResultsDict.items():
            ##if state does not change in given time: do something
            if state in Cond.keys():
                Cond[state].append([len(Nomin)/float(Denomin),len(Nomin)]) 
            else:
                Cond[state]=[len(Nomin)/float(Denomin),len(Nomin)]
        probabilities=[]
        for K, V in Cond.items():
            a=[K,V]
            probabilities.append([tuple(int(k) for k in K),V[0]])
        with open(self.saveResults+'\\Condprob{}.csv'.format(self.currentstate), "w") as Kap:
            Kap=csv.writer(Kap)
            for K, V in Cond.items():
                a=[K,V]
                Kap.writerow(a)
        self.probabilities=probabilities
        
    
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
    model = """
    model test
    compartment C1;
    C1 = 1e-12;
    species S1, S2;

    S1 = 9;
    S2 = 1;
    J1: 2S1+S2 -> 3S1; k1*S1^2*S2;
    J2: 3S1 -> 2S1+S2; k2*S1;
    k1 = 1;
    k2 = 6;
    end
    """
    mcsteps=200
    endTime=0.08
    fr=Firstreaction(model,mcsteps,endTime)
    fr.first_reaction()