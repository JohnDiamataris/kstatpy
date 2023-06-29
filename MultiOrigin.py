import numpy as np
import tellurium as te
import os
import pickle as pk
import csv
from roadrunner import Config
Config.setValue(Config.MAX_OUTPUT_ROWS, 1000000000)

class Multiorigin():
    def __init__(self):
        self.endTime=''
        self.stepRange=''
        self.i_states=''
        #self.allstates''=
        self.model=''
        self.results=''
        self.results2=''

    def multi_origin(self):
        try:
            z = te.loada(self.model)
        except:
            z = te.loadSBMLModel(self.model)
        z.integrator = 'gillespie'
        z.integrator.setValue('nonnegative', True)
        z.integrator.setValue('variable_step_size', True)
        #z.integrator.seed = 2
        s = z.simulate(0, self.endTime)
        i_state=list(s[0][1:])
        allstates=[i_state]
        i_states={}
        if not self.endTime == s[-1][0]:
            Config.setValue(Config.MAX_OUTPUT_ROWS, 10000000000)
            self.multi_origin()
        alldtval=self.endTime/len(s)        
        deci=1
        while alldtval*deci//1<= 0:
            deci *= 10
        deci = len(str(deci))
        alldtval=np.round(alldtval, deci)
        Tstep=alldtval
        for  st in range(1,self.stepRange):
            step=Tstep*st  
            timesteps=list(np.arange(0,s[-1][0]+step,step))
            s_index=1
            oldtime=0
            for k in timesteps:
                if k >= s.T[0][s_index]:
                    b=list(s.T[0]).index(s.T[0][s_index])
                    j_state=list(s[b][1:])
                    if not (j_state) in allstates:
                        allstates.append(j_state)
                    val=k-oldtime
                    decimal=1
                    while val*decimal//1<= 0:
                        decimal *= 10
                    decimal = len(str(decimal))+2
                    DTau=np.round(val, decimal)
                    i_state=tuple(i_state)
                    j_state=tuple(j_state)
                    if i_state in i_states.keys():
                        if j_state in i_states[i_state].keys():
                            if DTau in i_states[i_state][j_state].keys():
                                i_states[i_state][j_state][DTau][0]=(i_states[i_state][j_state][DTau][0]+1)
                            else:
                                i_states[i_state][j_state][DTau]=[1]
                        else:
                            i_states[i_state][j_state]={DTau:[1]}
                    else:
                        i_states[i_state]={j_state:{DTau:[1]}}

                    i_state=j_state
                    while s.T[0][s_index] < k and s_index <= len(timesteps) and k < s.T[0][-1]:
                        s_index += 1    
                    oldtime=k
        
        self.i_states= i_states
        #self.allstates= allstates
        self.statistic()

    def statistic(self):
    
        i_states=self.i_states
        #StatesList=self.allstates
        StatDict={}
        StatDict2={}
        
        for stat in i_states.items():
            stateSumOfj=0
            for i in stat[1].items():
                if not stat[0] in StatDict.keys():
                    StatDict[stat[0]]={}
                if not stat[0] in StatDict2.keys():
                    StatDict2[stat[0]]={}
                
                sumlenOfj=0
                lenOfJ=0
                for j in i[1].items():
                    lenOfJ+=j[1][0]
                sumlenOfj+=lenOfJ
                for j in i[1].items():
                    if not i[0] in StatDict[stat[0]].keys():
                        StatDict[stat[0]][i[0]]={}
                    if not i[0] in StatDict2[stat[0]].keys():
                        StatDict2[stat[0]][i[0]]=sumlenOfj
                        stateSumOfj+=sumlenOfj
                    if not j[0] in StatDict[stat[0]][i[0]].keys():
                        StatDict[stat[0]][i[0]][j[0]]=[j[1][0]/float(lenOfJ)]
        for i in StatDict2.items():
            all=sum(i[1].values())
            for p in i[1]:
                i[1][p]=i[1][p]/all 
            

        self.results = StatDict
        self.results2 = StatDict2
        resultspath=os.path.dirname(__file__)+'\\Results'
        with open(resultspath+'\\MultiOrigin_Stats.csv', "w") as Kap:
            Kap=csv.writer(Kap)
            for K, V in StatDict.items():
                a=[K,V]
                Kap.writerow(a)
        
          
    

if __name__ == "__main__":
    mo=Multiorigin()
    mo.multi_origin()
    