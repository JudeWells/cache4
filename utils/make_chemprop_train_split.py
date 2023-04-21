import pandas as pd
from utils.scaffold_split import load_train_test_data
train, test = load_train_test_data()

train = train[['smiles', 'target']]
test = test[['smiles', 'target']]
train.to_csv('ligand_files/chemprop_train.csv', index=False)
test.to_csv('ligand_files/chemprop_test.csv', index=False)