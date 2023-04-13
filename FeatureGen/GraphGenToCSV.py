# In[]
from pathlib import Path
from GraphGen import GraphInfoGen, get_label
import pandas as pd
import os
from pathlib import Path

#%%

filePath=Path()
print(filePath)
#%%
dataPath=filePath/'../ligand_files/CBLB_inhibitors.csv'
print(os.path.exists(dataPath))

df=pd.read_csv(dataPath)
df.info()
# %%
# GraphDF=pd.DataFrame(columns=['node features','bonds','bond features','target value'])
# %%
df=df.drop(columns=['sdf_file_ix','COMPOUND_NAMES','Reference'])
df.info()
df.head
# %%
df['Target']=df['IC50_range_nM'].apply(get_label)
df.head()
# %%
df['node features','bonds','bond features']=df['smiles'].apply(GraphInfoGen)
df.head()

# # %%
# for row_smile, row_target in zip(df['smiles'],df['IC50_range_nM']):
#     GraphRow=[]
#     nodeFeatures,bond,bond_f=GraphInfoGen(row_smile)
#     GraphRow.append(nodeFeatures)
#     GraphRow.append(bond)
#     GraphRow.append(bond_f)
#     GraphRow.append(get_label(row_target))
#     print(GraphRow)
#     GraphDF.loc[len(df.index)]=GraphRow
# # %%
# GraphDF.info()
# # %%
# GraphDF.to_csv('GraphDataSet.csv')
# # df['node features','bonds','bond features']=df.map)
# # %%

# %%
df.to_csv('GraphDataSet.csv')
# %%
df.info()
# %%
