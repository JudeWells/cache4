""" 

This note book contains functions which will generate the graph from a smiles string 


"""
# In[]
import os, re

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
import numpy as np

import torch
import torch_geometric.data as data
from torch_geometric.utils.convert import to_networkx

import networkx as nx

import time

from atomFeatureGen import Atom_features
from bondFeatureGen import bond_Feature

testSMILES='C[C@@H](c1cccc(c1)c1c2ccccc2[nH]n1)Sc1nncn1C'

def get_label(d):
    if d=="101-1000":
        return(550.5)
    if d=="1001-5000":
        return(3000.5)
    if d==">5000":
        return(5000)
    if d=="<100":
        return(100)
    if d=="101-300":
        return(200.5)
    if d=="301-1000":
        return(650.5)
    if d=="1001-3000":
        return(2000.5)
    if d=="3001-10000":
        return(6500.5)
    if d==">10000":
        return(10000) 

def MolGen(SMILES):
    #This function will take a Smiles string as input and will output a mol structure
    # start=time.time()
    mol=Chem.MolFromSmiles(SMILES)
    mol=Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    # AllChem.MMFFOptimizeMolecule(mol)
    AllChem.ComputeGasteigerCharges(mol)
    # end=time.time()-start
    # print(end)
    return mol

def nodeFeatureList(mol):
    #this takes as input a mol file and outputs a list of all of the nodes
    node_feature_list=[]
    for atom in mol.GetAtoms():
        node_feature_list.append(Atom_features(atom))
    return node_feature_list

def bondFeatures(mol):
    #this takes as input a mol file and outputs a list of the sorce and destination nodes and a list of bond features
    sourceAtomList=[]
    destAtomList=[]
    bondFeaturesList=[]
    for bond in mol.GetBonds():
        u=bond.GetBeginAtomIdx()
        v=bond.GetEndAtomIdx()
        bondFeature=bond_Feature(bond)
        sourceAtomList.extend([u,v])
        destAtomList.extend([v,u])
        bondFeaturesList.extend([bondFeature,bondFeature])
    atomSourceDestList=[sourceAtomList,destAtomList]
    return atomSourceDestList,bondFeaturesList

def GraphInfoGen(SMILES):
    #this function generates the relevent infomation to creater a pytorch Geometric
    mol=MolGen(SMILES)
    node_features=nodeFeatureList(mol)
    bonds,bond_f=bondFeatures(mol)
    node_features=torch.tensor(node_features,dtype=torch.float32)
    bonds=torch.tensor(bonds,dtype=torch.int32)
    bond_f=torch.tensor(bond_f,dtype=torch.float32)
    
    return node_features,bonds,bond_f

def GraphGen(SMILES,Target_Val):
    #This function generates a data object for pytorch Geometric 
    mol=MolGen(SMILES)
    node_features=nodeFeatureList(mol)
    bonds,bond_f=bondFeatures(mol)
    node_features=torch.tensor(node_features,dtype=torch.float32)
    bonds=torch.tensor(bonds,dtype=torch.int32)
    bond_f=torch.tensor(bond_f,dtype=torch.float32)

    return data.Data(x=node_features,edge_index=bonds,edge_attr=bond_f,y=Target_Val)

# %%
molTest=MolGen(testSMILES)
molTest
# Chem.MolToXYZ(molTest)
# testNodeList=nodeFeatureList(molTest)
# testBonds,testBond_f=bondFeatures(molTest)
# print(Chem.MolToMolBlock(molTest))
# # %%
# graph=GraphGen(testSMILES)
# print(graph)

# # %%
# vis = to_networkx(graph)

# # node_labels = graph.y.numpy()
# node_labels = graph.y
# import matplotlib.pyplot as plt
# plt.figure(1,figsize=(15,13)) 
# nx.draw(vis, cmap=plt.get_cmap('Set3'),node_size=70,linewidths=6)
# # nx.draw(vis, cmap=plt.get_cmap('Set3'),node_size=70,linewidths=6)
# plt.show()

# %%
