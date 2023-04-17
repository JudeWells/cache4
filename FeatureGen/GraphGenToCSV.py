# In[]
from pathlib import Path
from GraphGen import GraphInfoGen, get_label
import pandas as pd
import os
from pathlib import Path

filePath=Path()
print(filePath)
dataPath=filePath/'../ligand_files/CBLB_inhibitors.csv'
print(os.path.exists(dataPath))

df=pd.read_csv(dataPath)
df.info()

# %%
df=df.drop(columns=['sdf_file_ix','COMPOUND_NAMES','Reference'])
df.info()
df.head
# %%
df['Target']=df['IC50_range_nM'].apply(get_label)
df.head()
# %%
# df['smiles'].apply(GraphInfoGen)
# %%
df['tempstor']=df['smiles'].apply(GraphInfoGen,True)
df.info()

df.to_csv('GraphDataSet_additionalInfo.csv')
# %%
df.info()
# %%
