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
from process_amine_filter_infinisee import get_smiles_from_sdf_or_mol2, extract_molnames_from_sdf

def get_properties_from_sdf(sdf_path):
    """
    Extracts the following from the sdf file:
    Code (enamine ID)
    LOGP
    MW
    :param sdf_path:
    :return:
    """
    suppl = Chem.SDMolSupplier(sdf_path)
    mols = [m for m in suppl]
    names = [m.GetProp('Code') for m in mols]
    logP = [m.GetProp('LOGP') for m in mols]
    MW = [m.GetProp('MW') for m in mols]
    smiles = [Chem.MolToSmiles(m) for m in mols]
    return names, logP, MW, smiles

def main():
    sdf_path='final_selection/wim_best_poses_with_codes.sdf'
    en_ids, logPs, MWs, smiles = get_properties_from_sdf(sdf_path)
    new_rows = []
    for en_id, logP, MW, smile in zip(en_ids, logPs, MWs, smiles):
        new_row = {
            'Participant ID':'1972',
            'rank':'',
            'smiles':smile,
            'Enamine catalog ID':en_id,
            'mol_weight':MW,
            'logP':logP,
            'PAINS passed [y/n]':'',
            'Method':'wim pharmacophore',
            'MMGBSA': '',
            'Wim Recc.':1,
            'Jude Recc.':'',
            'Vina Score':'',
            'Pharmaco score':'',
        }
        new_rows.append(new_row)
    new_df = pd.DataFrame(new_rows)
    new_df.to_csv(sdf_path.replace('.sdf', '.csv'), index=False)


if __name__ == '__main__':
    main()