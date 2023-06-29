import os
from collections import OrderedDict
import matplotlib.pyplot as plt
import itertools as it
#from itertools import izip
from collections import Counter
import scipy
from scipy.integrate import simps
#import numpy as np
import pickle as pk
import csv
import re
#from lifelines import *
import scipy.optimize 
import scipy.stats as ss
from scipy.stats import expon
#import tellurium as te
import numpy as np
class MDAnalysis():
    
    def __init__(self):
        self.directorymain=''
        self.Timestep=1.
        self.systemdict={}
        self.allreactions=[]
        self.histoarrays=[]
        self.allstates=[]
    
    def main(self):
        #size=0
        system1=OrderedDict()
        allreactions=[]
        #time_int=[] ## list of time intervals
        allstates=[]
        for filename in self.directorymain[0]:
            with open(filename,'r')as f20:
                for line in f20:
                    if line.strip().startswith('#') == False:
                        a=line.strip() 
                        break             
            with open(filename,'r')as f:
                
                ZeroState=[]
                c = [-1 for x in range(len(re.split(';|,|\*|\n|\t|[\s]+',a)))]
                ZeroState.extend(c[:])
                steps = [ZeroState]
                # if (self.Timestep) == 0.:
                #     steps=[ZeroState]
                # else:
                #     steps=[ZeroState] ## first state to begin with
                flagstep=[] ## the step where an event occured
                events=[] ## each item in list is the original state and the state it went to (event).   
                timecount=-1
                for line in f:
                    
                    # size+=len(line)
                    # if size%1000 == 0:
                    #     size
                    if line.strip().startswith('#')==False:
                        if line.strip() == '':
                            line=next(f)
                            continue
                        steps.append(re.split(';|,|\*|\n|\t|[\s]+',line.strip()))
                        #steps.append(re.split(';|,|\*|\n|\t+',line.strip()))
                        if (self.Timestep) == 0.:
                            stepline0=[int(x) for x in steps[0][1:]]
                            stepline1=[int(x) for x in steps[1][1:]]
                            if stepline0 != stepline1:
                            ## when not equal then an event has occured and the calculations begin.
                                if stepline0 == ZeroState[1:]:
                                    flagstep = [int(float(x)) for x in steps[1]]
                                    del steps[0]
                                    continue
                                events.append([flagstep,steps[1]]) ## appends an event to the list event
                                time_old = float(flagstep[0]) ## time just before the event
                                time_new = float(steps[1][0]) ## time just after the event
                                
                                time_inter = time_new-time_old          
                                # start of the event
                                event_start = [int(x) for x in events[-1:][0][0][1:]]
            
                                # end of the event
                                event_end = [int(x) for x in steps[1][1:]]
                                flagstep = steps[1]

                                react= np.array(event_end)-np.array(event_start) ## reaction vector of current event
                                reaction = react.tolist()
                                if not reaction in allreactions:
                                    allreactions.append(reaction) ## appends a vector in the list if it isn't already in the list. 
                                
                                    
                                if not event_start in allstates:
                                    allstates.append(event_start)
                                    ind_oldstate = len(allstates)-1
                                else:
                                    ind_oldstate = allstates.index(event_start)
                                if not event_end in allstates:
                                    allstates.append(event_end)
                                    ind_newstate = len(allstates)-1
                                else:
                                    ind_newstate = allstates.index(event_end)                
                
                                key1=tuple(event_start)
                                value1=[time_inter,ind_newstate]
                                if key1 in system1.keys():  
                                    system1[key1].append(value1)  
                                else:  
                                    system1[key1] = [value1]
                                ## a dictionary with the states as keys, and values lists with the time spent in the state and the state it went to.
                            del steps[0]
                        else:
                            timecount+=1
                            stepline0=[int(x) for x in steps[0]]
                            stepline1=[int(x) for x in steps[1]]
                            if stepline0 != stepline1: ## when not equal then an event has occured and the calculations begin.
                                if stepline0 == ZeroState:
                                    flagstep = stepline1
                                    del steps[0]
                                    continue
                                events.append([flagstep,stepline1]) ## appends an event to the list event
                                time_inter = timecount*self.Timestep ## time it stayed in current state untill an event occured
                                ##time_int.append(time_inter)
                                timecount=0
                                # start of the event
                                event_start =  [int(x) for x in events[-1:][0][0]]
            
                                # end of the event
                                event_end = stepline1
                                flagstep = stepline1

                                react= np.array(event_end)-np.array(event_start) ## reaction vector of current event
                                reaction = react.tolist()
                                if not reaction in allreactions:
                                    allreactions.append(reaction) ## appends a vector in the list if it isn't already in the list. 
                                
                                    
                                if not event_start in allstates:
                                    allstates.append(event_start)
                                    ind_oldstate = len(allstates)-1
                                else:
                                    ind_oldstate = allstates.index(event_start)
                                if not event_end in allstates:
                                    allstates.append(event_end)
                                    ind_newstate = len(allstates)-1
                                else:
                                    ind_newstate = allstates.index(event_end)                
                
                                key1=tuple(event_start)
                                value1=[time_inter,ind_newstate]
                                if key1 in system1.keys():  
                                    system1[key1].append(value1)  
                                else:  
                                    system1[key1] = [value1]
                                ## a dictionary with the states as keys, and values lists with the time spent in the state and the state it went to.
                            del steps[0]
            
                self.systemdict=system1
                self.allreactions=allreactions
                self.allstates=allstates
                # Results = os.path.dirname(os.path.abspath(__file__))+'\\Results\\'
                # if not os.path.exists(Results):
                #     os.mkdir(Results)
                # fullpath2=os.path.join(Results,'Allstates.csv')
                # with open (fullpath2,'w') as f23:
                #     wr=csv.writer(f23)
                #     wr.writerows(self.allstates)  
                print ("finished file {}".format(filename))
        del events
       # del allstates
            
        for k,v in self.systemdict.items():
            self.histoarrays.append([k,v])
                
if __name__ == "__main__":
    x=MDAnalysis()
    x.main()