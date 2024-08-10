import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def make_smiles_canonical(df, smiles_col):
    df['canonical_smiles'] = df[smiles_col].apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x)))
    return df

def get_fingerprints(df, RADIUS=3, NBITS=2048):
    return [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), RADIUS, nBits=NBITS) for smi in df["canonical_smiles"]]

def visualise_molecules_that_match_on_enamine_id(quotes_df, original_df):
    mols_to_draw = []
    for i, row in quotes_df.iterrows():
        enamine_id = row['Catalog ID']
        original_match = original_df[original_df['Enamine catalog ID'] == enamine_id]
        if len(original_match) == 0:
            continue
        else:
            match_mol = Chem.MolFromSmiles(original_match.iloc[0]['smiles'])
            mols_to_draw.append(Chem.MolFromSmiles(row['New_SMILES']))
            mols_to_draw.append(match_mol)
    draw_mols(mols_to_draw)

def visualise_molecules_that_match_on_compound_id(quotes_df, original_df):
    mols_to_draw = []
    for i, row in quotes_df.iterrows():
        rank = int(row['ID'].split('-')[-1])
        original_match = original_df[original_df['rank'] == rank]
        if len(original_match) == 0:
            continue
        else:
            assert len(original_match) ==1
            match_mol = Chem.MolFromSmiles(original_match.iloc[0]['smiles'])
            mols_to_draw.append(Chem.MolFromSmiles(row['New_SMILES']))
            mols_to_draw.append(match_mol)
    draw_mols(mols_to_draw)


def find_most_similar_molecule(quotes_df, original_df):
    """
    Take the canonical smiles string from the quotes df
    and find the index of the most similar molecule in the
    original df (also using canonical smiles)
    add the most similar molecule's index in a new column
    and also include the tanimoto similarity score
    """
    quotes_fingerprints = get_fingerprints(quotes_df)
    original_fingerprints = get_fingerprints(original_df)
    similarity_matrix = np.zeros((len(quotes_fingerprints), len(original_fingerprints)))
    for i in range(len(quotes_fingerprints)):
        for j in range(len(original_fingerprints)):
            sim = DataStructs.FingerprintSimilarity(quotes_fingerprints[i], original_fingerprints[j])
            similarity_matrix[i, j] = sim

    for i, row in quotes_df.iterrows():
        similarity_scores = similarity_matrix[i, :]
        quotes_df.loc[i, 'most_similar_index'] = similarity_scores.argmax()
        quotes_df.loc[i, 'tanimoto_similarity'] = similarity_scores.max()
    return quotes_df


def draw_mols(mols):
    img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(200, 200))
    img.show()

if __name__=="__main__":
    local_df = pd.read_csv('final_selection/FINAL_SELECTIONS_148_Compound_selection_sample_file_challenge_4.csv')
    original_df = pd.read_csv('final_selection/WEBSITE_DOWNLOAD_23_AUG_1972_challenge-4_compound_selection.csv')
    original_df = make_smiles_canonical(original_df, 'smiles')

    quotes_df = pd.read_csv('final_selection/quotes_cache_4_1972.csv')
    for i, row in quotes_df.iterrows():
        if row.isnull().Original_SMILES:
            quotes_df.loc[i, 'Original_SMILES'] = row['New_SMILES']
    quotes_df = make_smiles_canonical(quotes_df, 'Original_SMILES')
    quotes_df = find_most_similar_molecule(quotes_df, original_df)

    quotes_df['compound_rank'] = quotes_df['ID'].apply(lambda x: int(x.split('-')[-1]))
    quotes_df.set_index('compound_rank', inplace=True, drop=False)
    original_df.set_index('rank', inplace=True, drop=False)
    joined = quotes_df.join(original_df[['Enamine catalog ID']])
    joined.set_index('Enamine catalog ID', inplace=True, drop=False)
    local_df.set_index('Enamine catalog ID', inplace=True, drop=False)
    joined = joined.join(local_df, lsuffix='_q', rsuffix='_local')
    bp=1


    visualise_molecules_that_match_on_enamine_id(quotes_df, original_df)
    visualise_molecules_that_match_on_compound_id(quotes_df, original_df)
    quotes_df.set_index('canonical_smiles', inplace=True, drop=False)
    original_df.set_index('canonical_smiles', inplace=True, drop=False)
    print(f"len(quotes_df) = {len(quotes_df)}")
    print(f"len(original_df) = {len(original_df)}")
    print(f"intersection canonical smiles = {len(set(original_df.canonical_smiles).intersection(set(quotes_df.canonical_smiles)))}")
    print(f"intersection enamine ids = {len(set(original_df['Enamine catalog ID']).intersection(set(quotes_df['Catalog ID'])))}")
    bp=1
    # combine the columns joined on the original smiles string
