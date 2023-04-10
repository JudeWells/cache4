import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools

sdf_file = "ligand_files/CBLB_inhibitors_vsF.sdf"

# Read SDF file using RDKit's SDMolSupplier
suppl = Chem.SDMolSupplier(sdf_file)

# Create an empty DataFrame with desired columns
columns = ["sdf_file_ix", "COMPOUND_NAMES", "Reference", "IC50_range_nM", "smiles"]
df = pd.DataFrame(columns=columns)

# Iterate through the molecules in the SDF file
rows = []
for i, mol in enumerate(suppl):
    if mol is not None:
        data = {}
        data["sdf_file_ix"] = i

        # Extract properties
        for prop in columns[1:]:
            if mol.HasProp(prop):
                data[prop] = mol.GetProp(prop)

        # Append the extracted data to the DataFrame
        rows.append(data)
df = pd.DataFrame(rows)

# Save the DataFrame to a CSV file
df.to_csv("ligand_files/CBLB_inhibitors.csv", index=False)