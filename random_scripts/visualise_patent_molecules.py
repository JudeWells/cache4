'''
Created by Jude Wells 2023-07-29
Draw the patent molecules in 2D
'''
from rdkit import Chem
from rdkit.Chem import Draw

def main():
    smiles_path = 'patented_mols/patented_mols.smi'
    with open(smiles_path) as f:
        smiles = [l.strip() for l in f.readlines()]
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    img = Draw.MolsToGridImage(mols, molsPerRow=int(len(mols)**0.5), subImgSize=(200, 200))
    img.save('patented_mols/patented_mols.png')

if __name__=="__main__":
    main()