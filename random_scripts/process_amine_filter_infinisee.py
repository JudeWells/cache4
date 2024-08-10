'''
Created by Jude Wells 2023-07-31
loads the sdf file containing the infinisee
molecules that have been filtered by the amine
and extracts the following:
Enamine ID
MMGBSA dG Bind
MW
logP
also generates a smiles string for each molecule
'''

import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def get_smiles_from_sdf_or_mol2(f_path):
    if f_path.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(f_path)
        mols = [m for m in suppl]
    elif f_path.endswith('.mol2'):
        suppl = Chem.MolFromMol2File(f_path)
        mols = [suppl]
    else:
        raise ValueError('File must be sdf or mol2')
    return [Chem.MolToSmiles(m) for m in  mols]

def extract_molnames_from_sdf(sdf_path):
    suppl = Chem.SDMolSupplier(sdf_path)
    mols = [m for m in suppl]
    return [m.GetProp('_Name') for m in mols]


def main():
    sdf_path='final_selection/infinisee_top2500_amine_filter.sdf'
    mmgbsa_path='infinisee/mmgbsa_infinisee_real_8gcy_top_2p5k_results.csv'
    names = extract_molnames_from_sdf(sdf_path)
    smiles = get_smiles_from_sdf_or_mol2(sdf_path)
    name2smiles = dict(zip(names, smiles))
    df = pd.read_csv(mmgbsa_path)
    df = df.dropna(subset=['MMGBSA dG Bind'])
    new_rows = []
    for k,v in name2smiles.items():
        df_row = df[df['Title']==k]
        if df_row.shape[0] == 0:
            print('No match for {}'.format(k))
            continue
        mmgbsa_result = df_row.iloc[0]['MMGBSA dG Bind']
        logP = df_row.iloc[0]['BIOSOLVEIT.LOGP']
        MW = df_row.iloc[0]['BIOSOLVEIT.MOLECULAR_WEIGHT']
        new_row = {
            'Participant ID':'1972',
            'rank':'',
            'smiles':v,
            'Enamine catalog ID':k,
            'mol_weight':MW,
            'logP':logP,
            'PAINS passed [y/n]':'',
            'Method':'infinisee REAL + MMGBSA  + amine',
            'MMGBSA':mmgbsa_result,
            'Wim Recc.':'',
            'Jude Recc.':'',
            'Vina Score':'',
            'Pharmaco score':'',
        }
        new_rows.append(new_row)
    new_df = pd.DataFrame(new_rows)
    new_df.to_csv(sdf_path.replace('.sdf', '.csv'), index=False)


if __name__ == '__main__':
    main()