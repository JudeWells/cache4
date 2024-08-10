import pandas as pd

"""
Desired mix:
22 - wim pharmacophore
10 - fingerprints
30 - infinisee
40 - vina
10 - classifier
30 - infinisee patent1
"""



def drop_bads(df, method):
    method_df = df[df.Method==method]
    assert len(method_df) > 1
    drop_indices = method_df[df['bad_apples'].isin(['High pScore'])]
    df.drop(drop_indices.index, inplace=True)
    method_df = df[df.Method == method]
    drop_indices = method_df[method_df.logP > 4.5]
    df.drop(drop_indices.index, inplace=True)
    return df


def keep_only_best_method_df(df, method, n):
    method_df = df[df.Method==method]
    method_df.sort_values(by='MMGBSA', inplace=True)
    mmgbsa_cut = method_df.iloc[n]['MMGBSA']
    drop_indices = method_df[method_df['MMGBSA'] > mmgbsa_cut]
    df.drop(drop_indices.index, inplace=True)
    return df

def main():
    df = pd.read_csv('final_selection/Compound_selection_sample_file_challenge_4.csv')
    methods = ['infinisee REAL + MMGBSA + amine',
                   'classifier REAL + MMGBSA + amine',
                    'vina + fingerprints',
                    'vina + MMGBSA',
               'infinisee patent1 + MMGBSA'
               ]
    keep_n = [31, 10, 10, 40, 31]
    for method in methods:
        df = drop_bads(df, method)

    for method, keep_n in zip(methods, keep_n):
        if method == 'vina + fingerprints':
            continue
        df = keep_only_best_method_df(df, method, keep_n)
    df.to_csv('final_selection/Compound_selection_sample_file_challenge_4_filtered.csv', index=False)



if __name__=="__main__":
    main()
