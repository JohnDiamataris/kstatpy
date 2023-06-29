import pickle as pk
import os
import numpy as np
import csv 
from decimal import Decimal
from VizualizeReactions import Reactions
from statsmodels.stats.weightstats import DescrStatsW

class WeightedAverage():
    def __init__(self):
        self.KappaAndError=''
        self.allreactions=''
        self.VizReactions={}
        self.modelspecies=[]
        

    def viz_reactions(self):
        vz=Reactions()
        vz.KappaAndError=self.KappaAndError
        vz.allreactions=self.allreactions
        vz.modelspecies=self.modelspecies
        vz.reactions()
        self.VizReactions=vz.VizReactions
        pass
    
    def weightedAv(self,Kaydict,StochKaydict):
        ''' calculates and returns the value of k as the weighted average and it's error '''
        filepath=os.path.dirname(os.path.abspath(__file__))+"\\Results"
        AllDict=[Kaydict,StochKaydict]
        KayandErr={}  ###---Kparam, Error, Cparam, error---######
        
        for KayDict in AllDict:
            for K in KayDict:
                """ allarray=np.array(KayDict[K])
                sumalln=sum(1/(allarray[:,1]**2))
                sumalln1=sum(allarray[:,2])
                weights=1/(allarray[:,1]**2)
                weights1=allarray[:,2]/sumalln1
                k_list=allarray[:,0]
                n_eff=sum(weights)**2/sum(weights**2)
                kk= sum([i*j for i,j in zip(k_list,weights)])/sumalln
                kk1=sum([i*j for i,j in zip(k_list,weights1)]) """
                Wi=[]
                Kapas=[]
                # calculation of Sum(1/error**2)
                for W in KayDict[K]:
                    Wi.append(1/W[1]**2)
                    Kapas.append(W[0]/W[1]**2)
                    # if np.isnan(W[1]):
                    #     Wi.append(1/W[1]**2)#####################################
                    # else: Wi.append(1/W[1]**2)
                sumallWi=np.nansum(Wi)
                # error of weighted average is sqrt(Sum(1/error**2))
                # calculating the weighted average of the k's gives us the approximated rate constant
                #Kapas=[]
                #n=0
                #for I in KayDict[K]:
                    # if not np.all(np.isnan(I[:-1])):
                    #     if np.isnan(I[1]):
                    #         #Wi.append(1/I[1]**2)##############################
                    #         Kapa += (((1/I[1]**2)*I[0])/sumallWi)#############
                    #     else:
                    #         #Wi.append(1/I[1]**2)
                    #Kapas.append(I[0]/I[1]**2)
                    #n+=I[2]
                #kkk=np.average(k_list,weights=weights)
                #weighted_stats = DescrStatsW(k_list, weights=weights, ddof=len(k_list))
                kapa = sum(Kapas)/sumallWi
                SE_k=np.sqrt(1/sumallWi)
                
                logVark=(1/kapa**2)*SE_k**2
                #CI=sorted([kapa+1.96*SE_k,kapa-1.96*SE_k])
                logCI=[np.log(kapa)+1.96*np.sqrt(logVark),np.log(kapa)-1.96*np.sqrt(logVark)]
                CI=sorted([np.exp(i) for i in logCI])
                if K in KayandErr.keys():
                    KayandErr[K].append(['{:.3e}'.format(kapa),'{:.3e}'.format(CI[0])+' - '+'{:.3e}'.format(CI[1])])
                else:
                    KayandErr[K]=[['{:.3e}'.format(kapa),'{:.3e}'.format(CI[0])+' - '+'{:.3e}'.format(CI[1])]]
        
              
        self.KappaAndError=KayandErr
        self.viz_reactions()
        with open(filepath+'\\WeightedAverage.csv','a')as fopen2:
            fopen2=csv.writer(fopen2)
            fopen2.writerow(['reaction','k  CI - c  CI'])
            fopen2.writerows(KayandErr.items())
           