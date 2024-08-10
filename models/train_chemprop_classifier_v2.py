import subprocess
import json


def run_chemprop_command(command):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {' '.join(command)}")
        print(result.stderr)
        return None
    return result.stdout


# Training command
train_command = [
    "chemprop", "train",
    "--data-path", "ligand_files/actives_decoys_drugbank_discodivers_ligands_drop_bad.csv",
    "-t", "classification",
    "--save-dir", "chemprop_2024_aug_test_new_dataset",
    "--split-type", "scaffold_balanced",
    "--target-columns", "target",
    "--ensemble-size", "1",
    "--num-folds", "3",
    "--metrics", "accuracy", "f1", "binary-mcc",
    "--split-sizes", "0.80", "0.1", "0.1",
    "--smiles-columns", "smiles",
    "--show-individual-scores"
]

if __name__ == '__main__':
    print("Starting Chemprop training...")
    train_output = run_chemprop_command(train_command)

    if train_output:
        print("Training completed. Output:")
        print(train_output)

        # Parse the JSON output to get mean and std scores
        try:
            scores = json.loads(train_output.split('\n')[-2])  # Assuming the last line is empty
            mean_scores = scores['mean']
            std_scores = scores['std']

            print("\nMean Scores:")
            for metric, value in mean_scores.items():
                print(f"{metric}: {value}")

            print("\nStandard Deviation Scores:")
            for metric, value in std_scores.items():
                print(f"{metric}: {value}")
        except json.JSONDecodeError:
            print("Could not parse score output. Please check the full output above.")