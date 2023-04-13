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
df['tempstor']=df['smiles'].apply(GraphInfoGen)
df.info()
# %%
print(len(df['tempstor'][1]))


# %%
df['node features']=df['tempstor']
# %%
df['bonds']=df['tempstor'][1]
df['bond features']=df['tempstor'][2]

# %%
df.to_csv('GraphDataSet.csv')
# %%
df.info()
# %%
