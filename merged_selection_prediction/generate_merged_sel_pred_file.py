import pandas as pd

df = pd.read_csv("merged_selection_prediction/231208_merged_selections_CACHE4.csv")
mmgbsa = pd.read_csv("merged_selection_prediction/mmgbsa_results_merged_selections_cache4_dg_bind_only.csv")
mmgbsa = mmgbsa[['MMGBSA dG Bind', 'Merged_ID']].dropna().set_index('Merged_ID', drop=False)
joined = df.join(mmgbsa, on='Merged_ID', how='left', lsuffix='_orig', rsuffix='_mmgbsa')
joined.sort_values('MMGBSA dG Bind', inplace=True)
chemprop_pred = pd.read_csv("merged_selection_prediction/chemprop_predictions_231208_merged_selections_CACHE4.csv")
chemprop_pred.set_index('Merged_ID', inplace=True, drop=False)
joined = joined.join(chemprop_pred, on='Merged_ID_orig', how='left', rsuffix='_chemprop')
joined = joined.rename(columns={"pred_0": "chemprop_pred"})
print(joined[['MMGBSA dG Bind', 'chemprop_pred']].corr())
joined.to_csv("merged_selection_prediction/merged_sel_pred_mmgbsa_chemprop.csv", index=False)

