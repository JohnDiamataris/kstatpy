import os
from collections import OrderedDict
import matplotlib.pyplot as plt
from collections import Counter
import scipy
from scipy.integrate import simps
import numpy as np
import pickle as pk
import csv
import re
import scipy.optimize 
import math
import scipy.stats as ss
from scipy.stats import expon
import tellurium as te

class SingleState():

    def __init__(self):
        self.directorymain=''
        self.Timestep=1.
        self.savehistoarray=''
        self.systemdict={}
        self.allreactions=[]
        self.histoarrays=[]
        self.allstates=[]
        self.probabilities=''
        self.singstate=0

  
    def single_state(self):
        savehistoarray = os.path.dirname(os.path.abspath(__file__))+'\\histoarrays'
        if not os.path.exists(savehistoarray):
            os.mkdir(savehistoarray)
        self.savehistoarray=savehistoarray
        cond_probability={}
        ltime_intervals=[]
        onlyTimeInt=[]
        events=[]   
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
                if (self.Timestep) == 0.:
                    steps=[ZeroState]
                else:
                    ZeroState.insert(0,-1)
                    steps=[ZeroState] ## first state to begin with
                flagstep=[]
                acc = -1
                for line in f:
                    acc += 1
                    if line in ['\n', '\r\n','']:
                        acc = 0
                        continue
                        
                    if self.Timestep != 0:
                        line = str(self.Timestep*acc)+','+line
                    if len(steps)==0:
                        steps.append(ZeroState)
                    steps.append(line.rstrip().split(','))
                    stplist=[int(x) for x in steps[0][1:]]
                          
                    if stplist != [int(x) for x in steps[1][1:]]:
                        if stplist == [int(x) for x in ZeroState[1:]]:
                            flagstep = steps[1]
                            del steps[0]
                            continue
                        events.append([flagstep,steps[1]])
                        if events[0][1] in [['\n'], ['\r\n']]:
                            del steps[:]
                            del events[:]
                            continue
                        time_old = float(flagstep[0])
                        time_new = float(steps[1][0])
                        time_inter = time_new-time_old
                        key1=tuple(events[0][1][1:])
                        key0=[int(x) for x in flagstep[1:]]
                        
                        # start of the event
                        event_start = [int(x) for x in events[-1:][0][0][1:]]
                        # end of the event
                        event_end = [int(x) for x in steps[1][1:]]
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
                        ltime_intervals.append([time_inter,ind_newstate])
                        onlyTimeInt.append(time_inter)
                        reaction=list(np.array(event_end)-np.array(event_start))
                        if not reaction in self.allreactions:
                            self.allreactions.append(reaction)
                        if key1 in cond_probability.keys():  
                            cond_probability[key1].append(time_inter)  
                        else:  
                            cond_probability[key1] = [time_inter]
                        del steps[0]
                        del events[:]
                        #pass
                        acc = 0
                        for line in f:
                            if line in ['\n', '\r\n']:
                                break
                    del steps[0]
                    del events[:]
                
                    

        with open('allstates.txt', "wb") as foallstates:
            pk.dump(allstates, foallstates)
        self.probabilities=[]
        for i,j in cond_probability.items():
            if len(i)==0:
                continue
            self.probabilities.append([i,len(j)/len(ltime_intervals)])
        self.histoarrays = [[tuple(key0),ltime_intervals]]   
        self.allstates = allstates    
        histoarray=os.path.join(self.savehistoarray, "ltimeintervals.txt")
        with open(histoarray, "wb") as f2:
            pk.dump([key0,ltime_intervals,self.probabilities], f2)
                
        with open('ltimeintervals.csv', 'w') as fo:
            datafl=csv.writer(fo)
            datafl.writerow(onlyTimeInt)
        with open('probabilities.csv', 'w') as f2:
            dataf2=csv.writer(f2)
            dataf2.writerow(self.probabilities)
        
   