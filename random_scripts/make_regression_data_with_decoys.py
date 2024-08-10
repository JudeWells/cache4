import pandas as pd
from random_scripts.plot_decoy_hit_properties import add_mol_properties
from rdkit import Chem
from rdkit.Chem import Draw

def visualise_molecules(df):
    # show the 2d molecular strcture of
    # a sample of 100 molecules
    # visualise the 2d structure using rdkit
    df = df[df.num_rings>5]
    df = df.sample(100)
    mols = [Chem.MolFromSmiles(x) for x in df.smiles]
    img = Draw.MolsToGridImage(mols, molsPerRow=10, subImgSize=(200, 200))
    img.save("ligand_files/decoy_molecules.png")
    img.show()



def main():
    df = pd.read_csv("ligand_files/all_decoys.csv")
    df.rename(columns={'SMILES': 'smiles'}, inplace=True)
    df = add_mol_properties(df)
    visualise_molecules(df)
    bp = 1

if __name__=="__main__":
    main()