import os
import matplotlib.pyplot as plt
from MeanTpopup import Ui_DialogMT
import pickle as pk
import scipy
from scipy.integrate import simps
from scipy.stats import expon
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines import ExponentialFitter
from BinwidthPopup import Ui_Dialog
from PyQt5 import QtCore, QtGui, QtWidgets
from lifelines.datasets import load_waltons
import csv

class MeanTime():

    def __init__(self):
        self.propbinwidth=''
        self.histoarray=[]
        self.binwidth=1E-12
        self.lengths=''
        self.lengthsandstate=''
        self.histoarrays=[]
        self.slider=''
        self.maxlen=''
        self.MeanTimeResults=''
        self.cellV=0
        self.allstates=''

    def MTpopup(self):
        Dialog = QtWidgets.QDialog()
        ui = Ui_DialogMT()
        ui.setupUi(Dialog)
        Dialog.exec()
        self.cellV=float(ui.lineEdit.text())
        
        #StochCorrec=ui.checkBox.isChecked()
        statelist=[]
        statelist.append(ui.lamda.isChecked())
        statelist.append(ui.trapezoid.isChecked())
        statelist.append(ui.km.isChecked())
        statelist.append(ui.arithm.isChecked())
        for ind,button in enumerate(statelist):
            if button is True:
                self.meanlife=ind

    def Binpopup(self):
        Dialog = QtWidgets.QDialog()
        ui = Ui_Dialog()
        ui.setupUi(Dialog)
        ui.proposedvalue.setText("{}".format(self.propbinwidth))
        Dialog.exec_()
        BinWidthValue=ui.lineEdit.text()
        if BinWidthValue == 0 or BinWidthValue == '':
            value = self.propbinwidth
        else:
            value  = BinWidthValue
        self.binwidth=float(value)

    def proposedbinwidth(self):
        histoarray1=[]
        for state in self.histoarrays:
            if len(state[1])>len(histoarray1):
                histoarray1=state[1][:]
        histlen=[]
        for i in histoarray1:
            histlen.append(i[0])
        self.propbinwidth=(max(histlen) - min(histlen))/(len(histlen))**0.5
        self.maxlen=len(histoarray1)
       
    
    def histogram(self):
        '''histogram creation and plotting '''
        
        histogramms=os.path.dirname(os.path.abspath(__file__))+'\\Results\\histogramms'
        if not os.path.exists(histogramms):
            os.mkdir(histogramms)
        
        MeanTimeResultsTxt = os.path.dirname(os.path.abspath(__file__))+'\\Results\\MeanTimeResultsTxt'
        if not os.path.exists(MeanTimeResultsTxt):
            os.mkdir(MeanTimeResultsTxt)
       
        KMcurves=os.path.dirname(os.path.abspath(__file__))+'\\Results\\KMcurves'
        if not os.path.exists(KMcurves):
            os.mkdir(KMcurves)
        
        lengths=[]
        lengthsandstate=[]
        self.MTpopup()
        MeanTimeResults=[]
        self.proposedbinwidth()
        self.Binpopup()
        
        for st in self.histoarrays:
            with open(histogramms +'//timelines{}.txt'.format(st[0]),'w')as hopen2:
                hopen2=csv.writer(hopen2)
                hopen2.writerow([st[0],st[1]])

            

         
        for state in self.histoarrays:
            nominator=[] ##list of states it goes to
            histoarray=[] ## list of lifetimes
            histoarray1=[] ##list of lifetimes and states it goes to
            formatd=[] ## state number
            histoarray2=[]  
            histoarray2.extend(state)
            histoarray1.extend(histoarray2[1])
            formatd.append(histoarray2[0])
            for i in histoarray1:
                histoarray.append(i[0])
                nominator.append(i[1])
            with open(histogramms +'//timelines{}.txt'.format(formatd),'w')as hopen2:
                hopen2=csv.writer(hopen2)
                hopen2.writerow(histoarray)
            self.histoarray=histoarray
            lengths.append(len(histoarray)) ## number of events
            lengthsandstate.append([formatd[0],len(histoarray)])
            denominator=float(len(nominator))
            trapezoid=0
            simpson=0   
            if len(histoarray)<2:
                continue
                 
            # n = plt.hist(histoarray, bins = np.arange(0, max(histoarray) + self.binwidth, self.binwidth)) ## number of events in each bin, bin start value
            # plt.close()
            # if len(histoarray)>(self.maxlen*(self.slider/100)):
            #     middle=[(x + y)*0.5 for x, y in zip(n[1][1:], n[1])] ## middle value of each bin for integration
            #     tpt=n[0]*middle
            #     mean_value_denom=scipy.integrate.simps(n[0], x=middle, dx=1, axis=-1, even='avg')
            #     mean_value_denom1=np.trapz(n[0], x=middle, axis=-1)
            #     mean_value_nom=scipy.integrate.simps(tpt, x=middle, dx=1, axis=-1, even='avg')
            #     mean_value_nom1=np.trapz(tpt, x=middle, axis=-1)
            #     trapezoid=mean_value_nom1/ mean_value_denom1
            #     simpson= mean_value_nom/ mean_value_denom
            # else:
            #     continue
            if len(histoarray)>(self.maxlen*(self.slider/100)):      
                '''Hazard analysis - finding the mean and variance of the Kaplan Meier function'''
                
                T=self.histoarray
                T0=[num for num in T if num]
                pol=min(T0)
                T1=[i/pol for i in T0]
                exf = ExponentialFitter()
                kmf = KaplanMeierFitter()
                kmf.fit(T)
                #exf.fit(T1,initial_point=[1e-24,])
                #exf.fit(T1)
                #exf._scipy_fit_options={'ftol': 1e-30, 'gtol': 1e-38}
                #exf._fit(initial_point=[1e-23,])
                #L=exf.lambda_
                #Lsum=exf.summary
                #LVar=exf.variance_matrix_._values[0][0]
                #L_Var=(L-Lsum._values[0][2])/1.96
                #L_Var=(L_Var**2)
                S=kmf.survival_function_.KM_estimate.tolist()
                S1=kmf.timeline
                # mean of K-M
                MeanArea=[]
                for m,q in enumerate(S[:-1]):
                    MeanArea.append((S1[m+1]-S1[m])*q)
                kmfMean=sum(MeanArea)
                state.append(kmfMean)
                #exf._bounds=([1.e-14,None])
                #exf._MIN_PARAMETER_VALUE=1.e-14
                exf.fit(T1)
                L=exf.lambda_*pol
                Lsum=exf.summary*pol
                LVar=exf.variance_matrix_._values[0][0]*pol**2
                
                L_Var=(L-Lsum._values[0][2])/1.96
                L_Var=(L_Var**2)
                Variance=0
                all=sum(kmf.event_table.T.values[1])
                for m,q in enumerate(S[:-1]):
                    ## events occured at time t(i)#########################
                    dead= kmf.event_table.T.values[1][m]
                    ## events that will occure at time t(i) and later#######
                    alive=kmf.event_table.T.values[-1][m] 
                    Variance += (sum(MeanArea[m:])**2 * dead)/float((alive*(alive-dead)))
                Variance=Variance * all/(all-1)
                Variances=[LVar,Variance]##  Lamda, KM
                average_time=float(sum(histoarray)/len(histoarray))
                #averages=[kmfMean]
                averages=['simpson',L,kmfMean,average_time]
                aaa =[[self.allstates[m],nominator.count(m)/denominator,nominator.count(m)] for m in set(nominator)] #state it goes to, probability and number of molecules.
                average_probability = [formatd[0],averages,aaa,Variances,denominator]  
                MeanTimeResults.append(average_probability)
                fullpath2=os.path.join(MeanTimeResultsTxt,'meantimeresults{}.txt'.format(formatd[0]))
                with open (fullpath2,'w') as f23:
                    f23.write('state,avereges,probabity to new state, Variances,observ')
                    f23.writelines( str(average_probability) )
                     
                        
        self.MeanTimeResults=MeanTimeResults
        self.lengths=lengths
        self.lengthsandstate=lengthsandstate
                    
            
                        

if __name__ == "__main__":
    a=MeanTime()
    a.histogram()