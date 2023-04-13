import pandas as pd
from rdkit import Chem
import json

path = 'ligand_files/CBLB_inhibitors.csv'
df = pd.read_csv(path)
dupes = df.duplicated(subset=['smiles']) 
df['Dupes'] = pd.Series(dupes)


# convert to canonised smiles
def canonise(df, path):
    smiles = df['smiles']
    cans_smiles = [Chem.MolToSmiles(Chem.MolFromSmiles(smi),True) for smi in smiles]
    smiles = pd.DataFrame(cans_smiles)
    print(df.head())
    df.to_csv(path)

def dupes(dupes):   
    return sum(dupes)

def dupes_to_json(df):
    dupe_rows = df.loc[df['Dupes'] == True]
    dupe_rows.to_csv('duplicated_rows_inhibitors.csv')

    seen = {} 
    count = 0  
    for i in dupe_rows['smiles']:
        if i in seen:
            seen[i] =  seen[i] + 1
            continue
        seen[i] = 2
    print(seen)
    with open('duplicated_row_counts.json', 'w') as fp:
        json.dump(seen, fp, sort_keys=True, indent=4)

dupes_to_json(df)    

                
