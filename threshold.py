from lifelines import ExponentialFitter
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
from scipy.stats import expon
from scipy.integrate import simps
import os
import numpy as np
import pickle as pk
from tkinter import filedialog as fd
import csv
import sys

Thresh=os.path.dirname(os.path.abspath(__file__))+'\\Thresh'
if not os.path.exists(Thresh):
    os.mkdir(Thresh)
#path = os.path.dirname(os.path.abspath(__file__))+'\\histoarrays' 
#if not os.path.exists(path):
    #os.mkdir(path)       
filenames = fd.askopenfilenames()
states=[] ## list of obs
threshold=[0.0,0.15,0.30,0.45,0.60,0.75]  ####   input 0 - 1
totalList=[]

for st in filenames:
    try:
        data = pk.load(open(st,'rb'))
        state=str(data[0]).strip(')(')
        for i in data[1]:
            states.append(i[0])
    except:
        with open(st,'r') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            mylist = list(csv_reader)
            csv.field_size_limit(100000000)
            obs= float(mylist[0][-1].strip(',][ '))
            states.append(obs)  
nmax= max(states)
for th in threshold:
    total=0
    for n in states:
        if n > nmax*th:                
            total+=n
    totalList.append([th,total])
with open(Thresh+os.sep+'MMETellurium_Longrun.csv'.format(th),'w') as csvfile:
    csvfile=csv.writer(csvfile)
    csvfile.writerow(totalList)
