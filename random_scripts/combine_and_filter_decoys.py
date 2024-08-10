import os

import pandas as pd

def load_decoy_finder_decoys():
    decoy_path = 'ligand_files/all_ligands.csv'
    df = pd.read_csv(decoy_path)
    df = df[df.target==0][['smiles']]
    df = df.rename(columns={'smiles': 'SMILES'})
    return df

def load_dude_decoys():
    dfs = []
    decoy_dir = 'ligand_files/sampled_hits_smiles/dude-decoys'
    for decoy_batch in os.listdir(decoy_dir):
        if 'decoys' not in decoy_batch:
            continue
        decoy_inner_path = os.path.join(*[decoy_dir, decoy_batch, 'decoys'])
        for decoy_file in os.listdir(decoy_inner_path):
            if decoy_file.endswith('.picked'):
                decoy_filepath = os.path.join(*[decoy_inner_path, decoy_file])
                df = pd.read_csv(decoy_filepath, sep='\t', header=None)
                df = df.iloc[1:, :-1]
                df.columns = ['SMILES', 'ZINC_ID', 'ZINC_Protonation_ID']
                dfs.append(df)
    df = pd.concat(dfs)
    return df

def combine_all_decoys():
    dude_decoys = load_dude_decoys()
    decoy_finder_decoys = load_decoy_finder_decoys()
    df = pd.concat([dude_decoys, decoy_finder_decoys])
    df = df.drop_duplicates(subset=['SMILES'])
    return df


def main():
    decoy_df = combine_all_decoys()
    decoy_df.to_csv('ligand_files/all_decoys.csv', index=False)
    bp=1

if __name__=="__main__":
    main()