import chemprop

arguments = [
    '--data_path', 'ligand_files/chemprop_train_v2.csv',
    '--dataset_type', 'classification',
    '--save_dir', 'chemprop_new_classifier',
    '--split_type', 'scaffold_balanced',
    '--target_columns', 'target',
    '--class_balance',
    '--ensemble_size', '1',
    '--num_folds', '5',
    '--extra_metrics', 'accuracy', 'f1', 'mcc',
    '--split_sizes', '0.9', '0.1', '0.0',
    '--smiles_column', 'smiles',
]

if __name__ == '__main__':
    args = chemprop.args.TrainArgs().parse_args(arguments)
    mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)
    print(mean_score)
    print(std_score)
