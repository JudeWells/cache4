"""
Created by Alice Cheng 2023-04-13
Functions for showing the number of duplicates and producing 
a new json showing duplicated entries. Then produces a dataset 
showing all molecules where the smiles entries has been seen more 
than once.
"""

import pandas as pd
from rdkit import Chem
import json

path = 'ligand_files/CBLB_inhibitors.csv'
df = pd.read_csv(path)
dupes = df.duplicated(subset=['smiles']) 
# create new column to store the duplicates in the dataset.
df['Dupes'] = pd.Series(dupes)

# dupes_bool = df['Dupes'].tolist()

# convert to canonised smiles
def canonise(df, path):
    smiles = df['smiles']
    cans_smiles = [Chem.MolToSmiles(Chem.MolFromSmiles(smi),True) for smi in smiles]
    smiles = pd.DataFrame(cans_smiles)
    print(df.head())
    df.to_csv(path)

# returns the total number of entries that are duplicated, 
# e.g. if CAS123 appears twice, there will be 1 dupeicated entry
# def dupes(dupes):   
#     return sum(dupes)

# creates a dictionary of all duplicated entries.
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

# for a smile string, return true if it has been seen more than once
def dupes_to_csv(smiles):
    seen_smi = dupes_seen(df).keys() 
    # print(smiles)
    if smiles in seen_smi:
        return True
    else:
        return False

 # iteratively applies dupes_to_csv to all smiles string column values
df['Dupes'] = df['smiles'].apply(dupes_to_csv)

# only show smiles strings where it has been seen more than once in the dataset.
duplicated = df.loc[df['Dupes'] == True]
duplicated.to_csv('duplicated_rows.csv', index=False)
     
            

                
