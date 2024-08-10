import pandas as pd
import numpy as np
import os
import sys
import shutil
from rdkit import Chem
"""
Created by Jude Wells 2021-07-31
Reads the molecules from sdf file
matches the result to the MMGBSA results
renames the sdf file to include the MMGBSA result
"""

def get_smiles_from_sdf_or_mol2(f_path):
    if f_path.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(f_path)
        mols = [m for m in suppl]
    elif f_path.endswith('.mol2'):
        suppl = Chem.MolFromMol2File(f_path)
        mols = [suppl]
    else:
        raise ValueError('File must be sdf or mol2')
    assert len(mols) == 1
    return Chem.MolToSmiles(mols[0])

def main():
    df_path = 'infinisee/mmgbsa_infinisee_real_8gcy_top_2p5k_results.csv'
    sdf_dir = 'infinisee/docked_poses_mmgbsa_top_150_infinisee_real'
    output_dir = os.path.join(sdf_dir, 'renamed')
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(df_path)
    df = df.dropna(subset=['MMGBSA dG Bind'])
    top150 = [f.split('.')[0] for f in os.listdir(sdf_dir) if f.endswith('.sdf') or f.endswith('.mol2')]
    df = df[df['Title'].isin(top150)]
    df.sort_values('MMGBSA dG Bind', inplace=True)
    df['rank_index'] = [str(n).zfill(3) for n in np.arange(df.shape[0])]
    for f in os.listdir(sdf_dir):
        if f.endswith('.mol2') or f.endswith('.sdf'):
            f_path = os.path.join(sdf_dir, f)
            df_row = df[df['Title']==f.split('.')[0]]
            if df_row.shape[0] == 0:
                print('No match for {}'.format(f))
                continue
            try:
                smiles = get_smiles_from_sdf_or_mol2(f_path)
            except:
                print('Could not get smiles for {}'.format(f))
                smiles = ''
                continue
            df.loc[df_row.index, 'smiles'] = smiles
            mmgbsa_result = str(int(round(df_row.iloc[0]['MMGBSA dG Bind'] * -1, 0)))
            rank_index = df_row.iloc[0]['rank_index']
            new_f = f"{rank_index}_dG{mmgbsa_result}_{f}"
            new_f_path = os.path.join(output_dir, new_f)
            shutil.copy(f_path, new_f_path)
            print('Renamed {} to {}'.format(f_path, new_f_path))
    df = df[['rank_index', 'Title', 'MMGBSA dG Bind', 'smiles', 'Prime MMGBSA ligand efficiency', 'source file', 'docking score']]
    df.to_csv(os.path.join(output_dir, df_path.split('/')[-1]), index=False)




if __name__=="__main__":
    main()