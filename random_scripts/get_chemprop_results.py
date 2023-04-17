import pandas as pd
from sklearn.metrics import balanced_accuracy_score, f1_score, matthews_corrcoef

preds = pd.read_csv('chemprop_preds.csv')
preds['target'] = preds['target'].apply(lambda x: 1 if x > 0.5 else 0)
true = pd.read_csv('ligand_files/chemprop_test.csv')

print('Bal Acc:', balanced_accuracy_score(true['target'], preds['target']))
print('MCC:', matthews_corrcoef(true['target'], preds['target']))
print('F1:', f1_score(true['target'], preds['target']))