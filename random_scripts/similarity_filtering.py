import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit.Chem import Draw

"""
Created by Jude Wells 2023-08-01
filters the final candidate selections
to ensure:
1) no duplicates
2) avoid highly similar molecules where possible
"""

def load_all_experimental_data():
    '''
    Create list of RDKIT mol objects for
    every molecule for which has already
    been tested experimentally
    :return:
    '''
    dfs = []
    ic50 = pd.read_csv('ligand_files/CBLB_inhibitors.csv')
    ic50['source'] = 'IC50'
    ic50 = ic50[['smiles', 'source']]
    dfs.append(ic50)
    patent = pd.read_csv('patented_mols/patented_mols.smi',header=None)
    patent.columns = ['smiles']
    patent['source'] = 'patent'
    dfs.append(patent)
    # get smiles representation of the mol2 file for the crystal 8gcy
    smiles_8gcy = 'FC(F)(F)c1c2c(cc(c1)C[N@H+]3C[C@H](CCC3)C)C(=O)N(c4cc(ccc4)[C@]5(C6=NN=CN6C)C[C@H](C5)C)C2'
    crystal = pd.DataFrame({'smiles':[smiles_8gcy], 'source':['8gcy']})
    dfs.append(crystal)
    return pd.concat(dfs)

def smiles_match_experimental(df, experimental):
    for i, row in df.iterrows():
        if row['smiles'] in experimental['smiles'].tolist():
            df.loc[i, 'smiles_match_experimental'] = 1
        else:
            df.loc[i, 'smiles_match_experimental'] = 0
    return df

def make_similarity_matrix(fp1, fp2):
    similarity_matrix = np.zeros((len(fp1), len(fp2)))
    for i in range(len(fp1)):
        for j in range(len(fp2)):
            sim = DataStructs.FingerprintSimilarity(fp1[i], fp2[j])
            similarity_matrix[i, j] = sim
    return similarity_matrix

def main():
    experimental = load_all_experimental_data()
    df = pd.read_csv('final_selection/Compound_selection_sample_file_challenge_4_filtered.csv')
    df = smiles_match_experimental(df, experimental)
    df_mols = df['smiles'].apply(lambda x: Chem.MolFromSmiles(x)).values
    radius = 3
    nbits = 2048
    df_fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits) for mol in df_mols]
    experimental_mols = experimental['smiles'].apply(lambda x: Chem.MolFromSmiles(x)).values
    experimental_fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits) for mol in experimental_mols]
    similarity_matrix = make_similarity_matrix(df_fingerprints, experimental_fingerprints)
    plt.hist(similarity_matrix.flatten(), bins=100)
    plt.title('Similarity between final selection and experimental molecules')
    plt.show()
    self_similarity_matrix = make_similarity_matrix(df_fingerprints, df_fingerprints)
    for i in range(len(self_similarity_matrix)):
        self_similarity_matrix[i,i] = 0
    plt.hist(self_similarity_matrix.flatten(), bins=100)
    plt.title('Similarity between final selection and final selection')
    plt.show()
    plt.imshow(similarity_matrix, interpolation='none')
    plt.title('Similarity between final selection and experimental molecules')
    plt.show()
    plt.imshow(self_similarity_matrix, interpolation='none')
    plt.title('Similarity between final selection with self')
    plt.show()
    df.reset_index(inplace=True, drop=True)
    assert len(df) == len(similarity_matrix)
    for i, sim_row in enumerate(self_similarity_matrix):
        sim_index = np.where(sim_row > 0.8)[0]
        if len(sim_index) > 0:
            df.loc[i, 'similar_to'] = '|'.join(df.loc[sim_index, 'Enamine catalog ID'].tolist())
            df.loc[i, 'similarity_scores'] = '|'.join([str(round(x, 2)) for x in sim_row[sim_index]])
    # sim_index = np.where(self_similarity_matrix > 0.8)
    # # draw pairs of molecules that are similar
    # sim_mols = []
    # for i in range(len(sim_index[0])):
    #     sim_mols.append(df_mols[sim_index[0][i]])
    #     sim_mols.append(experimental_mols[sim_index[1][i]])
    # legend = self_similarity_matrix[sim_index]
    # legend = [str(round(x, 2)) for x in legend]
    # double_legend = []
    # for i in range(len(legend)):
    #     double_legend.append(legend[i])
    #     enamine_ids = f"{df.loc[sim_index[0][i], 'Enamine catalog ID']} {df.loc[sim_index[1][i], 'Enamine catalog ID']}"
    #     double_legend.append(enamine_ids)
    # batch_size = 10
    # for i in range(0, len(sim_mols), batch_size):
    #     try:
    #         img = Draw.MolsToGridImage(sim_mols[i:i+batch_size], molsPerRow=2, subImgSize=(300, 500), legends=legend[i:i+batch_size])
    #         img.save(f'similar_mols_{i}.png')
    #     except:
    #         pass
    df.to_csv('final_selection/Compound_selection_sample_file_challenge_4_filtered.csv', index=False)


if __name__=="__main__":
    main()