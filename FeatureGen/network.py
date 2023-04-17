#In[]
import os, re
import numpy as np

import torch
import torch.nn as nn
from torch.nn.parameter import Parameter
import torch.nn.functional as F

from torch_geometric.nn import NNConv, global_add_pool
import torch_geometric.data as data
from torch_geometric.utils.convert import to_networkx

import networkx as nx

# %%

class testNet(torch.nn.Module):
    def __init__(self,numNodeFeatures,numEdgeFeatures):
        super().__init__()
        conv1_net=nn.Sequential(nn.Linear(numEdgeFeatures,32),
                                nn.ReLU(),
                                nn.Linear(32,numNodeFeatures*32))
        conv2_net=nn.Sequential(nn.Linear(numEdgeFeatures,32),
                                nn.ReLU(),
                                nn.Linear(32,32*16))
        self.conv1=NNConv(numNodeFeatures,32,conv1_net)
        self.conv2=NNConv(32,16,conv2_net)
        self.fc_1=nn.Linear(16,32)
        self.out=nn.Linear(32,1)

    def forward(self,data):
        batch,x,edge_index,edge_attr=(data.batch,data.x,data.edge_index,data.edge_attr)
        x=F.relu(self.conv1(x,edge_index,edge_attr))
        x=F.relu(self.conv2(x,edge_index,edge_attr))
        x=global_add_pool(x,batch)
        x=F.relu(self.fc_1(x))
        output=self.out(x)
        return output