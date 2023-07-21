arguments = [
    '--test_path', 'chemprop_test.csv',
    '--preds_path', 'chemprop_preds.csv',
    '--checkpoint_dir', 'chemprop_model5'
]

if __name__ == '__main__':
  args = chemprop.args.PredictArgs().parse_args(arguments)
  preds = chemprop.train.make_predictions(args=args)
