from PyQt5 import QtCore, QtGui, QtWidgets
from MeanTpopup import Ui_DialogMT
import os
import pickle as pk
import numpy as np
import math
import csv
from MeanTime import MeanTime
from Average import WeightedAverage
from VizualizeReactions import Reactions

class LeastSquares():

    def __init__(self):
        self.saveResults=''
        self.cellV=0
        self.meanlife=''
        self.lengths=''
        self.lengthsandstate=''
        self.slider=''
        self.histoarrays=[]
        self.MeanTimeResults=''
        self.allstates=''
        self.Results=''
        self.allreactions=[]
        self.VizReactions={}
        self.singlestate=0
        self.binwidth=''
        self.modelspecies=[]

    # def popup(self):
    #     Dialog = QtWidgets.QDialog()
    #     ui = Ui_Dialog()
    #     ui.setupUi(Dialog)
    #     Dialog.exec()
    #     self.cellV=float(ui.lineEdit.text())
       
    #     statelist=[]
    #     statelist.append(ui.simpson.isChecked())
    #     statelist.append(ui.trapezoid.isChecked())
    #     statelist.append(ui.km.isChecked())
    #     statelist.append(ui.arithm.isChecked())
    #     for ind,button in enumerate(statelist):
    #         if button is True:
    #             self.meanlife=ind
    
    def viz_reactions(self):
        vz=Reactions()
        vz.KappaAndError=self.KappaAndError
        vz.allreactions=self.allreactions
        vz.reactions()
        self.VizReactions=vz.VizReactions

    def meantime(self):
        mt=MeanTime()
        mt.histoarrays=self.histoarrays
        mt.slider=self.slider
        mt.allstates=self.allstates
        mt.histogram()
        self.binwidth=mt.binwidth
        self.MeanTimeResults=mt.MeanTimeResults
        self.lengths=mt.lengths
        self.lengthsandstate=mt.lengthsandstate
        self.cellV=mt.cellV
        self.meanlife=mt.meanlife


    def weighted_average(self,Kaydict=None,StochKaydict=None):
        av=WeightedAverage()
        av.allreactions=self.allreactions
        av.modelspecies=self.modelspecies
        av.weightedAv(Kaydict,StochKaydict)
        self.Results=av.KappaAndError
        self.VizReactions=av.VizReactions

    def Ksearch(self):
        saveResults = os.path.dirname(os.path.abspath(__file__))+'\\Results'
        if not os.path.exists(saveResults):
            os.mkdir(saveResults)
        self.meantime()
        self.saveResults=saveResults
        #self.popup()
        Kays={}
        StochKays={}#reaction: c,errorGreen,errorlamda,k,errorGreen,errorlamda
        allstates1=self.allstates ##all states that exist from whole simulation analysis
        NAvog=6.022140857e+23
        
        # maximum = int(max(self.lengths))
        # threshold = int(maximum * float(self.slider/100))
        # good_quality=[]
        # for i in self.lengthsandstate:
        #     if i[1]>=(threshold):
        #         good_quality.append(i[0])
        Kays1={}#'reaction': 'kay,Var,num,Deltat'
        StochKays1={}#'reaction': 'cee,Var,num,Deltat'
        for fn in self.MeanTimeResults:
            alldata=[]   
            reactionvector=[] ## reaction vector from current state to new state 
            stochvector=[]
            wnumber=[]#list of probability number of events and product of molecules 
            alldata.extend(fn) ## data of each state ((state number), average time, [state it goes to, probability, number of reactions],VarianceOfDt, numberofEvents ) 
            if alldata[3][self.meanlife-2]==0:
                continue
            # if not alldata[0] in good_quality:
            #     continue
            state=list(alldata[0]) ## current state
            for i in alldata[2]:##list [newstate,probability,number]
                try:
                    reactionvect=np.array(allstates1[i[0]])-np.array(state) ## reaction vector from current state to new state 
                except:
                    reactionvect=np.array(i[0])-np.array(state)
                reactionvector.append(reactionvect.tolist())
                reactants=[]
                for item in reactionvect:
                    if item<0:
                        reactants.append(abs(item))
                    else:
                        reactants.append(0)
                wnumber.append([i[1],i[2],np.prod((np.array(state)/(self.cellV*NAvog))**np.array(reactants))]) ## determenistic product                  
                #stochNumMultiplier=(np.prod((1/(self.cellV*NAvog))**np.array(reactants)))
                stochreact=[]## Stochastic  Product
                for k,n in zip(reactionvect, state):
                    if k<0:
                        stochreact.append(math.factorial(abs(n))/(math.factorial(abs(k))*math.factorial(abs(n)-abs(k))))
                    else:
                        stochreact.append(1)
                Stprod=1
                for p in stochreact:
                    Stprod *= p
                stochvector.append([i[1],i[2],Stprod])## list [prob, number, product]
            DeltaT1=[] ## mean lifetime in state 
            if self.meanlife == 1:       ###  placeholder
                DeltaT1.append(alldata[1][0])
            elif self.meanlife == 0:     ###  Lmean
                DeltaT1.append(alldata[1][1])
                Variance=alldata[3][0]
            elif self.meanlife == 2:       ## kaplan meier mean 
                DeltaT1.append(alldata[1][2])
                Variance=alldata[3][1]
            elif self.meanlife == 3:     ###  average time
                DeltaT1.append(alldata[1][3])
            DeltaT=DeltaT1[0]## DeltaT = 1/sumRi
            #Variance=alldata[3][self.meanlife-1]#*(self.cellV*NAvog)

            Nall=alldata[4]
        ##----C Parameter ------
            for item, reaction in zip(stochvector,reactionvector):
                Var=1
                kay=item[0]/(DeltaT*item[2])
                if Variance == alldata[3][0]:
                    Var=np.sqrt((kay/np.sqrt(Nall))**2 + (kay/(item[1]))**2*((item[1]*(1-item[0])))) ##errorLamda
                elif Variance == alldata[3][1]:
                    Var=np.sqrt((kay/DeltaT)**2 * Variance + (kay/(item[1]))**2*((item[1]*(1-item[0])))) ##errorGreen
                if tuple(reaction) in StochKays.keys(): 
                    StochKays[tuple(reaction)].append([kay,Var,item[1],DeltaT])  
                else:  
                    StochKays[tuple(reaction)] = [[kay,Var,item[1],DeltaT]] 
        ##----K Parameter ------    
            for item,reaction in zip(wnumber,reactionvector):
                Var=1
                kay=item[0]/(DeltaT*self.cellV*NAvog*item[2])
                if Variance == alldata[3][0]:
                    Var=np.sqrt((kay/np.sqrt(Nall))**2 + (kay/(item[1]))**2*((item[1]*(1-item[0]))))
                elif Variance == alldata[3][1]:
                    Var=np.sqrt((kay/DeltaT)**2 * Variance + (kay/(item[1]))**2*((item[1]*(1-item[0]))))
                if tuple(reaction) in Kays.keys(): 
                    Kays[tuple(reaction)].append([kay,Var,item[1],DeltaT]) 
                else:  
                    Kays[tuple(reaction)] = [[kay,Var,item[1],DeltaT]] 
        #####################  works with MDAnalysis
        for reaction in Kays.keys():
                reactants=[]
                products=[]
                vizreact=[]
                unitsum=0
                for factor,sp in zip(reaction,self.modelspecies):
                    if factor < 0:
                        reactants.append('{}{}'.format(abs(factor),sp.strip('[]')))
                        unitsum += abs(factor)
                    elif factor > 0:
                        products.append('{}{}'.format(factor,sp.strip('[]')))
                vizreact.append('+'.join(reactants))
                vizreact.append('+'.join(products))
                vizual='->'.join(vizreact)
                Kays1[vizual] = Kays[reaction]
                StochKays1[vizual] = StochKays[reaction]
        # for reaction in Kays.keys():
        #     reactants=[]
        #     products=[]
        #     vizreact=[]
        #     unitsum=0
        #     modspecies=self.modelspecies.split(',')
        #     for factor,sp in zip(reaction,modspecies):
        #         if factor < 0:
        #             reactants.append('{}{}'.format(abs(factor),sp))
        #             unitsum += abs(factor)
        #         elif factor > 0:
        #             products.append('{}{}'.format(factor,sp))
        #     vizreact.append('+'.join(reactants))
        #     vizreact.append('+'.join(products))
        #     vizual='->'.join(vizreact)
        #     Kays1[vizual] = Kays[reaction]
        #     StochKays1[vizual] = StochKays[reaction]
    
        with open(self.saveResults+'\\Kappa.csv', "w") as Kap:
            Kap=csv.writer(Kap)
            for K, V in Kays1.items():
                a=K,V
                Kap.writerow(a)
            for K, V in StochKays1.items():
                a=K,V
                Kap.writerow(a)
            # with open(self.saveResults+'\\Kappa.db', "wb") as pkKap:
            #     pk.dump(KayDict,pkKap)
        if self.singlestate==0:       
            self.weighted_average(Kays1,StochKays1)
        else:
            KayandErr={}
            for i,j in Kays1.items():
                if not i in KayandErr.keys():
                    KayandErr[i]=[[j[0][0],j[0][2]]]
                else: KayandErr[i].append([j[0][0],j[0][2]])
            for i,j in StochKays1.items():
                if not i in KayandErr.keys():
                    KayandErr[i]=[j[0][0],j[0][2]]
                else: KayandErr[i].append([j[0][0],j[0][2]])
            self.Results=KayandErr
            self.KappaAndError=KayandErr
            self.viz_reactions()
            with open(os.path.dirname(os.path.abspath(__file__))+'\\WeightedAverage.csv','a')as fopen2:
                fopen2=csv.writer(fopen2)
                fopen2.writerows(KayandErr)
       
if __name__ == "__main__":
    res=LeastSquares()
    res.Ksearch()