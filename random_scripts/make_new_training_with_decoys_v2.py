import pandas as pd
from rdkit import Chem

def load_and_combine_decoys_and_hits():
    decoys = pd.read_csv('ligand_files/all_decoys.csv')
    decoys['target'] = 0
    decoys.rename(columns={'SMILES': 'smiles'}, inplace=True)
    hits = pd.read_csv('ligand_files/all_ligands.csv')
    hits = hits[hits.target == 1]
    weak_hits = hits[hits.IC50_range_nM.isin(['>5000', '>10000'])].copy()
    weak_hits['target'] = 0
    hits = hits[~hits.IC50_range_nM.isin(['>5000', '>10000'])]
    combined = pd.concat([hits, weak_hits, decoys])
    return combined

def standardize_smiles(df):
    df['smiles'] = df['smiles'].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x)))
    return df


def main():
    df = load_and_combine_decoys_and_hits()
    df = standardize_smiles(df)
    df.to_csv('ligand_files/chemprop_train_v2.csv', index=False)

if __name__=="__main__":
    main()