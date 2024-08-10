import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from random_scripts.create_hits_subset_for_decoys import cluster_molecules
from models.train_chemprop_regression import convert_ic50_range
from rdkit.Chem.Draw import rdMolDraw2D

# Get the drawing options
opts = rdMolDraw2D.MolDrawOptions()

# Change the legend font size
opts.legendFontSize = 16


"""
Created by Jude Wells 2023-07-25
Script generates a handful of clusters
based on strucuture of the 900 mols
and measures the IC50 of the cluster
"""

def plot_ic50(df):
    """
    Dist plot of log IC50 for
    each cluster
    :return:
    """
    fig, axs = plt.subplots(len(df.cluster.unique()), 1, figsize=(10, 20))
    for i, cluster in enumerate(df.cluster.unique()):
        axs[i].set_title(f"Cluster {cluster} with {len(df[df.cluster == cluster])} mols")
        axs[i].set_xlim(df.log_ic50.min(), df.log_ic50.max())
        sns.distplot(df[df.cluster == cluster]['log_ic50'], hist=True, ax=axs[i], bins=30)
    plt.subplots_adjust(hspace=1)
    # reduce margins at top and bottom of plot
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.show()

def visualise_clusters(df, max_vis=36):
    cluster_results = {}
    for cluster in df.cluster.unique():
        cluster_mean = df[df.cluster == cluster].log_ic50.mean()
        cluster_results[cluster] = cluster_mean
        mols = np.array([Chem.MolFromSmiles(smi) for smi in df[df.cluster == cluster].smiles.values])
        IC50_scores = np.array(df[df.cluster == cluster].IC50_range_nM.values)
        max_vis = min(max_vis, len(mols))
        vis_index = np.random.choice(range(len(mols)), max_vis, replace=False)
        vis_mols = mols[vis_index]
        vis_ic50 = IC50_scores[vis_index]
        # generate legends for each molecule
        legends = [f'C{cluster} mean={round(cluster_mean, 3)}\n \nIC50={ic50}' for mol, ic50 in zip(vis_mols, vis_ic50)]
        img = Draw.MolsToGridImage(vis_mols, molsPerRow=6, subImgSize=(200, 200), legends=legends)
        img.show()
        bp=1

def main():
    df = pd.read_csv("ligand_files/CBLB_inhibitors.csv")
    df['ic50'] = df.IC50_range_nM.apply(convert_ic50_range)
    df['log_ic50'] = df.ic50.apply(lambda x: np.log10(x))
    df = cluster_molecules(df, smiles_column='smiles', n_clusters=6)
    # plot_ic50(df)
    visualise_clusters(df, max_vis=36)
    bp=1


if __name__=="__main__":
    main()