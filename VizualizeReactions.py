# import pickle as pk
# import os


class Reactions():

    def __init__(self):
        self.VizReactions=[]
        self.KappaAndError={}
        self.allreactions=''
        self.modelspecies=[]

    def reactions(self):
        vizual=[]
        if not self.modelspecies==[]:
            species = self.modelspecies
        else:
            species=[]
            for i in range(len(self.allreactions[0])):
                species.append('S{}'.format(i+1))
        reactionDict={}
        unitsdict={0:'mol/(l*sec)',1:'1/sec',2:'l/(mol*sec)',3:'l^2/(mol^2*sec)'}
        
       
        for reaction in self.allreactions:
            reactants=[]
            products=[]
            vizreact=[]
            unitsum=0
            for factor,sp in zip(reaction,species):
                if factor < 0:
                    reactants.append('{}{}'.format(abs(factor),sp.strip('[]')))
                    unitsum += abs(factor)
                elif factor > 0:
                    products.append('{}{}'.format(factor,sp.strip('[]')))
            vizreact.append('+'.join(reactants))
            vizreact.append('+'.join(products))
            vizual='->'.join(vizreact)
            if unitsum in unitsdict.keys():
                unit=unitsdict[unitsum]
            else: unit= 'unspecified'
            if vizual in reactionDict.keys():
                reactionDict[vizual].append([vizual,unit])
            else:
                reactionDict[vizual]=[vizual,unit]
            
        self.VizReactions=reactionDict
            

if __name__ == "__main__":
    r=Reactions()
    r.reactions()
             

 