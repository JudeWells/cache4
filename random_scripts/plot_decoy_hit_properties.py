import pandas as pd
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def add_mol_properties(df):
    """
    Adds the following columns to the df:
    - molecular weight
    - num_heavy_atoms
    - num_rings
    - num_rotatable_bonds
    - num_h_bond_donors
    - num_h_bond_acceptors
    :param df:
    :return:
    """
    df["mol"] = df.smiles.apply(lambda x: Chem.MolFromSmiles(x))
    df["mol_weight"] = df.mol.apply(lambda x: Descriptors.MolWt(x))
    df["num_heavy_atoms"] = df.mol.apply(lambda x: rdMolDescriptors.CalcNumHeavyAtoms(x))
    df["num_rings"] = df.mol.apply(lambda x: rdMolDescriptors.CalcNumRings(x))
    df["num_rotatable_bonds"] = df.mol.apply(lambda x: rdMolDescriptors.CalcNumRotatableBonds(x))
    df["num_h_bond_donors"] = df.mol.apply(lambda x: rdMolDescriptors.CalcNumHBD(x))
    df["num_h_bond_acceptors"] = df.mol.apply(lambda x: rdMolDescriptors.CalcNumHBA(x))
    return df

def plot_distributions(df):
    """
    Plots the distributions (histogram density plot)
    of both the hits (target=1) and the decoys (target=0)
    of the following properties:
    - molecular weight
    - num_heavy_atoms
    - num_rings
    - num_rotatable_bonds
    - num_h_bond_donors
    - num_h_bond_acceptors
    :param df:
    :return:
    """
    colnames = ["mol_weight", "num_heavy_atoms", "num_rings", "num_rotatable_bonds",
                "num_h_bond_donors", "num_h_bond_acceptors"]
    for colname in colnames:
        plt.figure()
        sns.distplot(df[df.target==0][colname], hist=True, label="decoys", color="skyblue",bins=30)
        sns.distplot(df[df.target==1][colname], hist=True, label="hits", color="red", bins=30)
        # sns.histplot(data1, color="skyblue", label="Data 1", kde=True, stat='density')
        # sns.histplot(data2, color="red", label="Data 2", kde=True, stat='density')
        plt.legend()
        plt.title(colname)
        plt.show()


def main():
    df = pd.read_csv("ligand_files/all_ligands.csv")
    df = df[df.target==1]
    decoy_df = pd.read_csv("ligand_files/all_decoys.csv")
    decoy_df["target"] = 0
    decoy_df = decoy_df.rename(columns={'SMILES': 'smiles'})
    df = pd.concat([df, decoy_df])
    df = add_mol_properties(df)
    plot_distributions(df)



if __name__=="__main__":
    main()