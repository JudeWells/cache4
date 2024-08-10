"""
Created by Jude Wells 2023-07-25
check if the relationship between
IC50 can be explained by any simple
molecular properties
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from models.train_chemprop_regression import convert_ic50_range
from random_scripts.plot_decoy_hit_properties import add_mol_properties

def train_and_eval_basic_regression(df):
    """
    Train a basic regression model
    on the data and evaluate it
    :param df:
    :return:
    """
    model = LinearRegression()
    features_all = ['mol_weight', 'num_rings', 'num_rotatable_bonds',
                  'num_h_bond_donors', 'num_h_bond_acceptors']

    features_rings = ['num_rings']

    features_weight = ['mol_weight']

    for features in [features_all, features_rings, features_weight]:
        print(f'\nFeatures: {features}')
        model.fit(df[features], df.log_ic50)
        for k,v in zip(features, model.coef_):
            print(k, round(v, 3))
        r2_score = model.score(df[features], df.log_ic50)
        print('r2_score', round(r2_score, 3))

def plot_hexgrid_heatmap_rings_vs_ic50(df):
    """
    Plot a heatmap of the number of rings vs the log ic50
    :param df:
    :return:
    """
    fig, ax = plt.subplots()
    # Creating the hexbin plot
    hb = ax.hexbin(df.num_rings, df.log_ic50, gridsize=10, cmap='inferno')

    # Adding a colorbar
    cb = plt.colorbar(hb, ax=ax)
    cb.set_label('counts')

    plt.xlabel('num_rings')
    plt.ylabel('log_ic50')
    plt.title('Hexbin plot with counts in each bin')
    plt.show()


def main():
    df = pd.read_csv("ligand_files/CBLB_inhibitors.csv")
    df['ic50'] = df.IC50_range_nM.apply(convert_ic50_range)
    df['log_ic50'] = df.ic50.apply(lambda x: np.log10(x))
    add_mol_properties(df)
    plot_hexgrid_heatmap_rings_vs_ic50(df)
    train_and_eval_basic_regression(df)

if __name__ == '__main__':
    main()
