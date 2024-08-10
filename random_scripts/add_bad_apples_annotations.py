import os
import pandas as pd

def load_all_bad_apples():
    bad_app_dir = 'bad_apples'
    dfs = []
    for fname in os.listdir(bad_app_dir):
        df = pd.read_csv(os.path.join(bad_app_dir, fname))
        dfs.append(df)
    df = pd.concat(dfs)
    return df

def get_badapp_val(apples_match):
    if 'High pScore.' in apples_match.advisory.values:
        return 'High pScore'
    elif 'Moderate pScore.' in apples_match.advisory.values:
        return 'Moderate pScore'
    elif 'Low pScore.' in apples_match.advisory.values:
        return 'Low pScore'
    else:
        return 'No pScore'

def main():
    df_path = 'final_selection/Compound_selection_sample_file_challenge_4.csv'
    apples = load_all_bad_apples()
    apples.rename(columns={c:c.strip() for c in apples.columns}, inplace=True)
    df = pd.read_csv('final_selection/Compound_selection_sample_file_challenge_4.csv')
    for i, row in df.iterrows():
        smiles = row['smiles']
        apples_match = apples[apples['mol_smi'] == smiles]
        if len(apples_match.mol_i.unique()) != 1:
            print(f'{len(apples_match.mol_i.unique())} for', smiles)
            bp=True
        bad_val = get_badapp_val(apples_match)
        df.loc[i, 'bad_apples'] = bad_val
    df[df['bad_apples'].isin(['High pScore'])].to_csv('final_selection/bad_apples_high_p_score_only.csv', index=False)
    df.to_csv(df_path, index=False)

if __name__=="__main__":
    main()