# In[]
import os, re
import numpy as np
import matplotlib.pyplot as plt

import torch
import torch.nn as nn
from torch.nn.parameter import Parameter
import torch.nn.functional as F

from torch_geometric.nn import NNConv, global_add_pool
import torch_geometric.data as data
from torch_geometric.utils.convert import to_networkx
from torch_geometric.loader import DataLoader

from torch.utils.data import random_split

import networkx as nx

from network import testNet
from GraphDataLoader import GraphDataset


device=torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)

# %%
dser=GraphDataset(ReactionFile='GraphDataSet.csv',PathToData='.')

#%%
trainSetSize=int(len(dser)*0.8)
valSetSize=int((len(dser)-trainSetSize)/2)
testSetSize=(len(dser)-trainSetSize)-valSetSize
print(len(dser))
trainData,valData,testData=random_split(dser,[trainSetSize,valSetSize,testSetSize])
print(len(trainData)+len(valData))

   
nodeFeatures,edgeFeatures=dser[1].num_features,dser[1].num_edge_features

# add data loader here
trainLoader=DataLoader(trainData,batch_size=32,shuffle=False)
valLoader=DataLoader(valData,batch_size=32,shuffle=False)
testLoader=DataLoader(testData,batch_size=32,shuffle=False)
# for batch in valLoader:
#     print(batch.y)
# %%
# network initilizaition 

net = testNet(nodeFeatures,edgeFeatures)

#optimixer 
optimizer=torch.optim.Adam(net.parameters(),lr=0.01)
net.to(device)

epochs=6

# Training  
# %%
for total_epochs in range(epochs):
    epoch_loss=0
    total_graphs=0
    net.train()
    for batch in trainLoader:
        batch.to(device)
        optimizer.zero_grad()
        output=net(batch)
        loss=F.mse_loss(output,batch.y[:,0])
        loss.backward()
        epoch_loss+=loss.item()
        total_graphs+=batch.num_graphs
        optimizer.step()
    train_avg_loss=epoch_loss/total_graphs
    val_loss=0
    total_graphs=0
    net.eval()
    for batch in valLoader:
        batch.to(device)
        output=net(batch)
        loss=F.mse_loss(output,batch.y[:,0])
        val_loss+=loss.item()
        total_graphs+=batch.num_graphs
    val_avg_loss=val_loss/total_graphs
    print(f"Epocs: {total_epochs} | epoch avg. loss: {train_avg_loss:.2f} | validation avg. loss: {val_avg_loss:.2f}")

# %%
net.eval()
predictions=[]
real=[]
for batch in testLoader:
    output=net(batch.to(device))
    predictions.append(output.detach().cpu().numpy())
    real.append(batch.y[:].detach().cpu().numpy())
real=np.concatenate(real)
predictions=np.concatenate(predictions)


# %%
plt.scatter(real[:],predictions[:])
plt.xlabel('actual IC_50_nM')
plt.ylabel('predicted IC_50_nM')
# plt.plot([0,100],[0,100])
# plt.xlim()
# %%
