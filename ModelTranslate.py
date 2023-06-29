#import libsbml
#from libsbml import *
import tellurium as te
#import numpy as np
import os
#te.setDefaultPlottingEngine('matplotlib')


class ModelManagement():
    def __init__(self):
        self.model=''
        self.savepath=os.path.dirname(__file__)+'\\ModelTranslate'
        if not os.path.exists(self.savepath):
            os.mkdir(self.savepath)

    
    def ExportSBML(self):#(model)
        '''Exports model as SBML, input model self.model and path to save to'''
        r = te.loada(self.model)
        sbmlmodel=r.getCurrentSBML()
        with open(self.savepath+'\\sbmlmodel.xml', 'w') as savefile:
            print (sbmlmodel,file=savefile)


    def ExportAntimony(self):#(model)
        '''Exports model as Antimony, input model name and path to save to'''
        r = te.loadSBMLModel(self.model)
        antimmodel=r.getCurrentAntimony()
        with open(self.savepath+'\\antimonymodel.ant','w') as savefile:
            print (antimmodel,file=savefile)



