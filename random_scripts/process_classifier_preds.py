'''
Created by Jude Wells 2023-07-31
loads the processed sdf file
into a format compatible with the final
submission
'''

import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors


def get_mols(sdf_path):
    suppl = Chem.SDMolSupplier(sdf_path)
    return [m for m in suppl]


def main():
    sdf_path='final_selection/classifier_739_amine_filter.sdf'
    mmgbsa_path='classifier_predictions/mmgbsa_classifier_preds_real.csv'
    mols = get_mols(sdf_path)
    names = [m.GetProp('_Name') for m in mols]
    smiles = [Chem.MolToSmiles(m) for m in mols]
    mol_weights = [Descriptors.MolWt(m) for m in mols]
    logps = [Descriptors.MolLogP(m) for m in mols]
    df = pd.read_csv(mmgbsa_path)
    df = df.dropna(subset=['MMGBSA dG Bind'])
    new_rows = []

    for name, smile, mol_weight, logp in zip(names, smiles, mol_weights, logps):
        one_mol = df[df['Title']==name]
        assert one_mol.shape[0] == 1
        enamine_id = one_mol['idnumber'].iloc[0]
        mmgbsa = one_mol['MMGBSA dG Bind'].iloc[0]
        new_row = {
            'Participant ID':'1972',
            'rank':'',
            'smiles':smile,
            'Enamine catalog ID':enamine_id,
            'mol_weight':mol_weight,
            'logP':logp,
            'PAINS passed [y/n]':'',
            'Method':'classifier REAL + MMGBSA  + amine',
            'MMGBSA':mmgbsa,
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