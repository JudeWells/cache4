from utils.scaffold_split import make_decoy_df, get_scaffold
import pandas as pd

def get_label(d):
    if d=="101-1000":
        return(550.5)
    if d=="1001-5000":
        return(3000.5)
    if d==">5000":
        return(5000)
    if d=="<100":
        return(100)
    if d=="101-300":
        return(200.5)
    if d=="301-1000":
        return(650.5)
    if d=="1001-3000":
        return(2000.5)
    if d=="3001-10000":
        return(6500.5)
    if d==">10000":
        return(10000)


def make_joined_dataset(hit_path="ligand_files/CBLB_inhibitors.csv",
                        decoy_path="ligand_files/decoys/found_decoys.sdf", save_path=None):
    decoys = make_decoy_df(decoy_path)
    decoys["ic_50_nM"] = 20000
    hits = pd.read_csv(hit_path)
    hits["ic_50_nM"] = hits["IC50_range_nM"].apply(get_label)
    hits["target"] = 1
    hits["scaffold"] = hits["smiles"].apply(get_scaffold)
    df = pd.concat([hits, decoys])
    if save_path is not None:
        df.to_csv(save_path, index=False)
    return df

bp=1