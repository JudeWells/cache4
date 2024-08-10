import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from models.train_chemprop_regression import convert_ic50_range
import seaborn as sns

def plot_glide_score_ic50(df, xcol='docking score'):
    ic50_bands = ['<100', '101-300', '101-1000', '301-1000', '1001-3000', '1001-5000', '3001-10000', '>5000', '>10000']
    fig, axs = plt.subplots(len(ic50_bands), 1, figsize=(10, 20))
    for i, band in enumerate(ic50_bands):
        axs[i].set_title(f"IC50 {band} nM")
        axs[i].set_xlim(df[xcol].min(), df[xcol].max())
        sns.distplot(df[df['IC50_range_nM'] == band][xcol], hist=True, ax=axs[i], bins=30)
    r2 = np.corrcoef(df[xcol], df['log_ic50'])[0, 1] ** 2
    plt.suptitle(f'{xcol} vs IC50 R2={round(r2, 3)}')
    plt.subplots_adjust(hspace=1)
    # reduce margins at top and bottom of plot
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.show()
    bp=1




def main():
    df_path = 'glide_docked_900_ic50/mmgbsa_900_ic50_results.csv'
    df = pd.read_csv(df_path)
    df = df[df['docking score'].notnull()]
    df['ic50'] = df.IC50_range_nM.apply(convert_ic50_range)
    df['log_ic50'] = df.ic50.apply(lambda x: np.log10(x))
    plot_glide_score_ic50(df, xcol='MMGBSA dG Bind')



if __name__=="__main__":
    main()