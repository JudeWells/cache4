import pandas as pd
from rdkit import Chem
import json

path = 'ligand_files/CBLB_inhibitors.csv'
df = pd.read_csv(path)
dupes = df.duplicated(subset=['smiles']) 
df['Dupes'] = pd.Series(dupes)

# dupes_bool = df['Dupes'].tolist()

# convert to canonised smiles
def canonise(df, path):
    smiles = df['smiles']
    cans_smiles = [Chem.MolToSmiles(Chem.MolFromSmiles(smi),True) for smi in smiles]
    smiles = pd.DataFrame(cans_smiles)
    print(df.head())
    df.to_csv(path)

# def dupes(dupes):   
#     return sum(dupes)

def dupes_seen(df):
    
    dupe_rows = df.loc[df['Dupes'] == True]
    seen = {} 
    for i in dupe_rows['smiles']:
        if i in seen:
            seen[i] =  seen[i] + 1
            continue
        seen[i] = 2
    # print(seen.keys())
    # with open('duplicated_row_counts.json', 'w') as fp:
    #     json.dump(seen, fp, sort_keys=True, indent=4)
    return seen

def dupes_to_csv(smiles):
    seen_smi = dupes_seen(df).keys() 
    # print(smiles)
    if smiles in seen_smi:
        return True
    else:
        return False

df['Dupes'] = df['smiles'].apply(dupes_to_csv)
duplicated = df.loc[df['Dupes'] == True]
duplicated.to_csv('duplicated_rows.csv', index=False)
     
            

                