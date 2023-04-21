import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import train_test_split


def get_scaffold(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
        return Chem.MolToSmiles(scaffold_generic)
    except:
        return None

def load_mol_objects_from_sdf(sdf_filepath):
    suppl = Chem.SDMolSupplier(sdf_filepath)
    return [mol for mol in suppl if mol is not None]

def make_decoy_df(sdf_filepath):
    decoys = load_mol_objects_from_sdf(sdf_filepath)
    decoy_df = pd.DataFrame()
    decoy_df["smiles"] = [Chem.MolToSmiles(mol) for mol in decoys]
    decoy_df["scaffold"] = decoy_df["smiles"].apply(get_scaffold)
    decoy_df["target"] = 0
    return decoy_df

def load_train_test_data(hits_path="ligand_files/CBLB_inhibitors.csv",
                         decoys_path="ligand_files/decoys/found_decoys.sdf", test_size=0.12, random_state=42):
    hits_df = pd.read_csv(hits_path)
    hits_df['target'] = 1
    hits_df["scaffold"] = hits_df["smiles"].apply(get_scaffold)
    decoy_df = make_decoy_df(decoys_path)
    full_df = pd.concat([hits_df, decoy_df])
    # Perform train-test split on unique scaffolds
    unique_scaffolds = full_df["scaffold"].unique()
    train_scaffolds, test_scaffolds = train_test_split(unique_scaffolds, test_size=test_size, random_state=random_state)
    # Split hit and decoy data using scaffolds
    train = full_df[full_df["scaffold"].isin(train_scaffolds)]
    test = full_df[full_df["scaffold"].isin(test_scaffolds)]
    print(f"train_hits: {train.shape}")
    print(f"test_hits: {test.shape}")
    print(f"train value counts: \n{train['target'].value_counts()}")
    print(f"test value counts: \n{test['target'].value_counts()}")
    print(f"proportion hits train: {train.target.mean()}")
    print(f"proportion hits test: {test.target.mean()}")
    return train, test


if __name__=="__main__":
    train, test = load_train_test_data()
    bp=1