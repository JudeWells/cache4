import os
import pandas as pd
import numpy as np
import chemprop

arguments = [
    '--dataset_type', 'regression',
    '--save_dir', 'chemprop_regression_logIC50_2_folds',
    '--split_type', 'scaffold_balanced',
    '--target_columns', 'log_ic50',
    '--num_folds' , '2',
    '--ensemble_size', '1',
    '--extra_metrics', 'r2', 'mae'
    '--split_sizes', '0.85', '0.15', '0.0',
]

def convert_ic50_range(ic50):
    if '-' in ic50:
        start, end = ic50.split('-')
        return (float(start) + float(end)) / 2
    elif ic50=='>5000':
        return 7500
    elif ic50=='<100':
        return 2
    elif ic50=='>10000':
        return 15000


def make_chemprop_regression_dataset(data_path):
    df = pd.read_csv('ligand_files/CBLB_inhibitors.csv')
    df['ic50'] = df.IC50_range_nM.apply(convert_ic50_range)
    df['log_ic50'] = df.ic50.apply(lambda x: np.log10(x))
    new_df = df[['smiles', 'log_ic50']]
    new_df = new_df.rename(columns={'log_ic50': 'target'})
    new_df.to_csv(data_path, index=False)


if __name__ == '__main__':
    data_path = 'ligand_files/chemprop_regression_no_split.csv'
    if not os.path.exists(data_path):
        make_chemprop_regression_dataset(data_path)
    arguments.append('--data_path')
    arguments.append(data_path)
    args = chemprop.args.TrainArgs().parse_args(arguments)
    mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)
    print(mean_score)
    print(std_score)
