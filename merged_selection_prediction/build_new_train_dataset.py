import csv
import os
import pandas as pd
from rdkit import Chem

def sdf_to_smiles_csv(sdf_path, csv_path):
    # Read the SDF file
    suppl = Chem.SDMolSupplier(sdf_path)

    # Open the CSV file for writing
    with open(csv_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)

        # Write the header
        csv_writer.writerow(['smiles', "inchikey", "mol_id", "sdf_file_index"])
        fail_counter = 0
        success_counter = 0
        # Iterate through the molecules in the SDF file
        for i, mol in enumerate(suppl):
            if mol is not None:
                # Generate SMILES string for the molecule
                smiles = Chem.MolToSmiles(mol)
                inchikey = Chem.MolToInchiKey(mol)
                catalog_id = mol.GetProp('Catalog ID') if mol.HasProp('Catalog ID') else ''
                # Write the SMILES string to the CSV file
                csv_writer.writerow([smiles, inchikey, catalog_id, i])
                success_counter += 1
            else:
                print("Warning: Skipped an invalid molecule")
                fail_counter += 1

    print(f"Conversion complete. CSV file saved at: {csv_path}")
    print(f"Successfully converted {success_counter} molecules")
    print(f"Failed to convert {fail_counter} molecules")

def combine_datasets():
    actives_and_decoys = pd.read_csv(actives_and_decoys_path)
    drugbank = pd.read_csv(drugbank_path)
    discovery_diversity = pd.read_csv(discovery_diversity_smiles_path)
    # Add a column to each dataset to indicate the source
    # if sdf_file_ix is present source is "known_actives"
    actives_and_decoys.loc[actives_and_decoys['sdf_file_ix'].notnull(), 'source'] = 'known_actives'
    actives_and_decoys.loc[actives_and_decoys['sdf_file_ix'].isnull(), 'source'] = 'decoy_finder'
    drugbank['source'] = 'drugbank'
    discovery_diversity['source'] = 'enamine_discovery_diversity'
    drugbank_reformatted = drugbank.rename(columns={
        'name': 'COMPOUND_NAMES',
    })[['smiles', 'COMPOUND_NAMES', 'source']]
    drugbank_reformatted['target'] = 0
    drugbank_reformatted['ic_50_nM'] = 20000

    discover_reformatted = discovery_diversity.rename(columns={
        "SMILES": "smiles",
        "mol_id": "COMPOUND_NAMES"
    })[['smiles', 'COMPOUND_NAMES', 'source']]
    discover_reformatted['target'] = 0
    discover_reformatted['ic_50_nM'] = 20000
    lit_actives = pd.read_csv(literature_path)
    combined = pd.concat([actives_and_decoys, drugbank_reformatted, discover_reformatted, lit_actives])
    combined.to_csv(combined_path, index=False)
    bp=1

    pass
def main():
    if not os.path.exists(discovery_diversity_smiles_path):
        sdf_to_smiles_csv(discovery_diversity_sdf, discovery_diversity_smiles_path)
    if not os.path.exists(combined_path):
        combine_datasets()







if __name__ == '__main__':
    actives_and_decoys_path = "ligand_files/all_ligands.csv"
    drugbank_path = "ligand_files/drugbank_smiles.csv"
    discovery_diversity_sdf = "ligand_files/Enamine_Discovery_Diversity_Set_50240_DDS-50_50240cmpds_20240729.sdf"
    discovery_diversity_smiles_path = "ligand_files/Enamine_Discovery_Diversity_Set_50240_DDS-50_50240cmpds_20240729.csv"
    patented_path = "patented_mols.smi"
    literature_path = "ligand_files/patented_mols/jw_found_literature_actives.csv"
    combined_path = "ligand_files/actives_decoys_drugbank_discodivers_ligands.csv"
    main()