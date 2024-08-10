'''
Created by Jude Wells 2023-07-29
Draw the patent molecules in 2D
'''
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

def main():
    df = pd.read_csv('ligand_files/CBLB_inhibitors.csv')[['IC50_range_nM','smiles']]
    df = df[df['IC50_range_nM']=='<100']
    smiles = df['smiles'].tolist()
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    legends = [f'{i} {smi[:10]}' for i, smi in enumerate(smiles)]
    img = Draw.MolsToGridImage(mols, molsPerRow=int(len(mols)**0.5), subImgSize=(300, 300), legends=legends)
    img.save('ligand_files/sub100_mols.png')

if __name__=="__main__":
    main()