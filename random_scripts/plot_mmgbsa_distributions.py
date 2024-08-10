import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
'''
Created by Jude Wells 2023-07-31
plots the distribution density of
MMGBSA scores for molecules generated
from different pipelines:
- infinisee on REAL
- chemprop classifier on REAL
- Vina docking on REAL
- Glide docking on HLL
'''

def main():
    path_dict = {
        'infinisee': 'infinisee/mmgbsa_infinisee_real_8gcy_top_2p5k_results.csv',
        'vina': 'vina/mmgbsa_vina_3500_results.csv',
        'chemprop': 'classifier_predictions/mmgbsa_classifier_preds_real.csv',
        'disco_divers_glide': 'schrodinger_output/mmgbsa_disco_divers_top_2k_results.csv'
    }
    # one subplot for each pipeline - each one with seaborn kdeplot
    fig, axs = plt.subplots(4, 1, figsize=(10, 26))
    colors = {1:'red', 2:'blue', 3:'green', 4:'purple'}
    for i, (pipeline, path) in enumerate(path_dict.items()):
        df = pd.read_csv(path)
        sns.distplot(df['MMGBSA dG Bind'], hist=True, ax=axs[1], bins=30, label=pipeline, color=colors[i+1])
        axs[1].legend()
        axs[i].set_title(pipeline)
    axs[1].legend()
    # set xlims to be the same for all subplots
    xlims = [ax.get_xlim() for ax in axs]
    xlims = [min([x[0] for x in xlims]), max([x[1] for x in xlims])]
    for ax in axs:
        ax.set_xlim(xlims)
    plt.show()

if __name__=="__main__":
    main()