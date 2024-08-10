import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

"""
Created by Jude Wells 2023-07-29
combine the predictions from the cluster
add some molecular propertties to the df
visualise the molecular properties of 
the top scorers
"""

def load_and_combine_preds():
    dir = 'classifier_predictions'
    dfs = []
    for subdir in os.listdir(dir):
        subdir_path = os.path.join(dir, subdir)
        if not os.path.isdir(subdir_path):
            continue
        for subsubdir in os.listdir(subdir_path):
            subsubdir_path = os.path.join(subdir_path, subsubdir)
            if not os.path.isdir(subsubdir_path):
                continue
            for f in os.listdir(subsubdir_path):
                if f.endswith('.csv'):
                    df_path = os.path.join(subsubdir_path, f)
                    one_df = pd.read_csv(df_path)
                    one_df['enamine_collection'] = subdir + '_REAL'
                    dfs.append(one_df)
    return pd.concat(dfs).reset_index(drop=True)

def draw_mols(mols):
    img = Draw.MolsToGridImage(mols, molsPerRow=int(len(mols)**0.5), subImgSize=(200, 200))
    img.show()

def draw_top_molecules(df, n, reverse=False):
    df = df.sort_values('pred', ascending=reverse).reset_index(drop=True)
    mols = [Chem.MolFromSmiles(smi) for smi in df.iloc[:n].smiles.values]
    draw_mols(mols)

def draw_random_molecules(df, n):
    df = df.sample(n)
    mols = [Chem.MolFromSmiles(smi) for smi in df.smiles.values]
    draw_mols(mols)

def main():
    df = load_and_combine_preds()
    create_sample_of_molecules(df, n_random=5000, n_top=7000)
    draw_random_molecules(df, 100)
    for reverse in [False, True]:
        draw_top_molecules(df, 100, reverse=reverse)

def create_sample_of_molecules(df, n_random, n_top):
    # create a sample of molecules
    top_sample_indices = df.sort_values('pred', ascending=False).iloc[:n_top].index
    random_sample_indices = df.sample(n_random).index
    sample_indices = np.concatenate([top_sample_indices, random_sample_indices])
    sampled_df = df.loc[sample_indices]
    sampled_df = sampled_df.sort_values('pred', ascending=False)
    sampled_df = sampled_df.drop_duplicates(subset=['smiles'])
    sampled_df = sampled_df[['smiles', 'idnumber']]
    sampled_df.to_csv('classifier_predictions/classifier_preds_real.csv', index=False)



if __name__=="__main__":
    main()