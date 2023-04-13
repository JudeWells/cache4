# %%
import networkx as nx
import numpy as np
import pandas as pd

import torch
from torch.utils.data import random_split ,Dataset #,DataLoader

import torch_geometric.data as data
# from torch_geometric.data import Dataset #this seems to brake it for some reason
from torch_geometric.loader import DataLoader
from torch_geometric.utils.convert import to_networkx

import networkx as nx


class GraphDataset(Dataset):
    def __init__(self,ReactionFile=None,PathToData='',target='Yield',Temp=None,Pressure=None,Time=None,AtomFeatures=None,BondFeatures=None):
        self.DFSMILES=pd.read_csv(str(PathToData)+'/'+str(ReactionFile))
        print(self.DFSMILES.info())
    def __len__(self):
        return len(self.DFSMILES)
    def __getitem__(self,idx):
        print(self.DFSMILES.loc[idx,'Target'])
        GraphList=eval(self.DFSMILES.loc[idx,'tempstor'])
        print(len(GraphList))
        combinedfeatures=torch.tensor(GraphList[0],dtype=torch.float32)
        bonds=torch.tensor(GraphList[1],dtype=torch.int64)
        bond_f=torch.tensor(GraphList[2],dtype=torch.float32)

        # combinedfeatures=torch.tensor(self.DFSMILES.loc[idx,'node features'],dtype=torch.float32)
        # bonds=torch.tensor(self.DFSMILES.loc[idx,'bonds'],dtype=torch.int64)
        # bond_f=torch.tensor(self.DFSMILES.loc[idx,'bond features'],dtype=torch.float32)
        label=self.DFSMILES.loc[idx,'Target']
        graph=data.Data(x=combinedfeatures,edge_index=bonds,edge_attr=bond_f,y=torch.tensor([[label]],dtype=torch.float32))
        
        return graph
    

# %%
dser=GraphDataset(ReactionFile='GraphDataSet.csv',PathToData='.')

# %%
trainSetSize=int(len(dser)*0.001)
testSetSize=len(dser)-trainSetSize
print(len(dser))
trainData,testData=random_split(dser,[trainSetSize,testSetSize])
print(len(trainData)+len(testData))
print(dser[30])
print(dser[30].num_features)
print(dser[30].num_edge_features)

# %%
graph=dser[30]
# %%
vis = to_networkx(graph)

# node_labels = graph.y.numpy()
node_labels = graph.y
import matplotlib.pyplot as plt
plt.figure(1,figsize=(15,13)) 
nx.draw(vis, cmap=plt.get_cmap('Set3'),node_size=70,linewidths=6)
# nx.draw(vis, cmap=plt.get_cmap('Set3'),node_size=70,linewidths=6)
plt.show()

# %%
