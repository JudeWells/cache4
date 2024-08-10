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
from random_scripts.process_classifier_preds import get_mols

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
    sdf_path='infinisee/infinisee_patent1/infinisee_patent_top_mmgbsa.sdf'
    mmgbsa_path='infinisee/infinisee_patent1/mmgbsa_results_infinisee_patent1_all_results.csv'
    df = pd.read_csv(mmgbsa_path)
    df = df.dropna(subset=['MMGBSA dG Bind'])
    new_rows = []
    mols = get_mols(sdf_path)
    names = [m.GetProp('_Name') for m in mols]
    smiles = [Chem.MolToSmiles(m) for m in mols]
    mol_weights = [Descriptors.MolWt(m) for m in mols]
    logps = [Descriptors.MolLogP(m) for m in mols]


    for name, smile, mol_weight, logp in zip(names, smiles, mol_weights, logps):
        one_mol = df[df['Title']==name]
        assert one_mol.shape[0] == 1
        enamine_id = one_mol['Title'].iloc[0]
        mmgbsa = one_mol['MMGBSA dG Bind'].iloc[0]
        new_row = {
            'Participant ID':'1972',
            'rank':'',
            'smiles':smile,
            'Enamine catalog ID':enamine_id,
            'mol_weight':mol_weight,
            'logP':logp,
            'PAINS passed [y/n]':'',
            'Method':'infinisee patent1 + MMGBSA',
            'MMGBSA':mmgbsa,
            'Wim Recc.':'',
            'Jude Recc.':'',
            'Vina Score':'',
            'Pharmaco score':'',
        }
        new_rows.append(new_row)
    new_df = pd.DataFrame(new_rows)
    new_df.sort_values(by='MMGBSA', inplace=True)
    new_df = new_df.iloc[:50]
    new_df.to_csv('final_selection/infinisee_patent_1.csv', index=False)


if __name__ == '__main__':
    main()