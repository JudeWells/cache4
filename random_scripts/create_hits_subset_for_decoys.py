"""
Created by Jude Wells 2023-07-22
Purpose of script is to select
a representative sample 45 of the hit
ligands so that we can create decoys 
using DUD-E
"""
import os

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.DataStructs import DiceSimilarity
from sklearn.cluster import AgglomerativeClustering
import pandas as pd
import numpy as np


def calculate_similarity_matrix(fp_list):
    size = len(fp_list)
    similarity_matrix = np.zeros((size, size))

    for i in range(size):
        for j in range(i + 1, size):
            similarity = DiceSimilarity(fp_list[i], fp_list[j])
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity

    return similarity_matrix


def smiles_to_fps(smiles_list, n_bits=1024):
    fps = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
        fps.append(fp)
    return fps


def cluster_molecules(df, smiles_column, n_clusters):
    smiles = df[smiles_column].tolist()
    fps = smiles_to_fps(smiles)
    similarity_matrix = calculate_similarity_matrix(fps)
    clustering = AgglomerativeClustering(n_clusters=n_clusters, affinity='precomputed', linkage='average')
    clusters = clustering.fit_predict(1 - similarity_matrix)
    df['cluster'] = clusters
    return df

def generate_smiles_files(df, n_files=3, n_per_file=45, savedir='ligand_files/sampled_hits_smiles'):
    """
    Samples one molecule from each cluster and saves a file with the smiles of the sampled molecules
    """
    os.makedirs(savedir, exist_ok=True)
    smiles_lists = []
    for cluster in df["cluster"].unique():
        cluster_df = df[df["cluster"]==cluster]
        smiles_lists.append(cluster_df.sample(1)["smiles"].tolist()[0])
    for i in range(n_files):
        with open(os.path.join(savedir, f"ic50_sampled_smiles_{i}.smi"), 'w') as f:
            f.write("\n".join(smiles_lists[i*n_per_file:(i+1)*n_per_file]))


def main():
    np.random.seed(13)
    df = pd.read_csv("ligand_files/CBLB_inhibitors.csv")
    n_per_file = 45
    n_files = 4
    df = cluster_molecules(df, "smiles", n_clusters=n_files*n_per_file)
    generate_smiles_files(df, n_files=n_files, n_per_file=n_per_file, savedir='ligand_files/sampled_hits_smiles_v2')
    bp=1


if __name__=="__main__":
    main()