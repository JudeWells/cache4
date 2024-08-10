import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import Descriptors

def get_mol_props(one_row):
    pdb_fname = one_row["Title"].replace('.sdf', '.pdb')
    pdb_path = os.path.join('vina/candidates_5000', pdb_fname)
    # load molecule from pdb file using rdkit
    mol = Chem.MolFromPDBFile(pdb_path)
    smiles = Chem.MolToSmiles(mol)
    mol_weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    mmgbsa = one_row['MMGBSA dG Bind']
    return smiles, mol_weight, logp, mmgbsa

def add_vina_fingerprints_best(df):
    new_preds = pd.read_csv('vina/fingerprints_best.txt', sep=' ', header=None)
    new_preds.columns = ['Title', 'fingerprint_score']
    new_preds['MMGBSA dG Bind'] = -999
    df = pd.concat([df, new_preds])
    return df


def main():
    df = pd.read_csv('vina/mmgbsa_vina_3500_results.csv')
    df = add_vina_fingerprints_best(df)
    lookup = pd.read_csv('vina/lookup_candidates_5000.csv')
    df.dropna(subset=['MMGBSA dG Bind'], inplace=True)
    df.sort_values(by='MMGBSA dG Bind', inplace=True)
    df.reset_index(inplace=True, drop=True)
    selected = df.iloc[:200]
    new_rows = []
    for i, row in selected.iterrows():
        row = row.to_dict()
        try:
            smiles, mol_weight, logp, mmgbsa = get_mol_props(row)
        except:
            print('failed to get mol props for', row['Title'])
            continue
        lookup_match = lookup[lookup['file'] == row['Title'].replace('.sdf', '.pdb')]
        assert lookup_match.shape[0] == 1
        vina_score = lookup_match['docking_score'].iloc[0]
        enamine_id = lookup_match['enamine_real_id'].iloc[0]
        new_row = {
            'Participant ID': '1972',
            'rank': '',
            'smiles': smiles,
            'Enamine catalog ID': enamine_id,
            'mol_weight': mol_weight,
            'logP': logp,
            'PAINS passed [y/n]': '',
            'Method': 'vina + MMGBSA',
            'MMGBSA': mmgbsa,
            'Wim Recc.': '',
            'Jude Recc.': '',
            'Vina Score': vina_score,
            'Pharmaco score': '',
        }
        if mmgbsa == -999:
            new_row['Method'] = 'vina + fingerprints'
            new_row['Wim Recc.'] = 1
        new_rows.append(new_row)
    new_df = pd.DataFrame(new_rows)
    new_df.to_csv('final_selection/vina_mmgbsa_selected.csv', index=False)
    bp=True

if __name__=="__main__":
    main()