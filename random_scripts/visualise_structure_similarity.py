import pandas as pd
import numpy as np
from rdkit.Chem import AllChem
from rdkit import Chem, DataStructs
import matplotlib.pyplot as plt
import seaborn as sns

"""
Created by Jude Wells 2023-07-20
Create a similarity matrix for all ligands
where similarity is defined as tanimoto similarity
between the morgan fingerprints of the ligands
"""

def get_fingerprints(df, RADIUS=3, NBITS=2048):
    return [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), RADIUS, nBits=NBITS) for smi in df["smiles"]]

def plot_fingerprints(fingerprints):
    # reduce margins at top and bottom of plot
    plt.figure(figsize=(5, 20))
    plt.imshow(fingerprints, interpolation='none')
    plt.show()

def create_similarity_matrix(fingerprints):
    similarity_matrix = np.zeros((len(fingerprints), len(fingerprints)))
    for i in range(len(fingerprints)):
        for j in range(i+1):
            sim = DataStructs.FingerprintSimilarity(fingerprints[i], fingerprints[j])
            similarity_matrix[i, j] = sim
            similarity_matrix[j, i] = sim
    return similarity_matrix

def main():
    df_path = "ligand_files/all_ligands.csv"
    df = pd.read_csv(df_path)
    fingerprints = get_fingerprints(df)
    similarity_mat = create_similarity_matrix(fingerprints)
    plt.figure(figsize=(10, 10))
    sns.heatmap(similarity_mat, cmap='viridis')
    plt.title('Tanimoto Similarity Matrix')
    plt.show()
    bp=1


if __name__=="__main__":
    main()