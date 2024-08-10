'''
Created by Jude Wells 2023-07-29
load the infiisee matched molecules
that were matched to the patent molecule
and then remove any that have already
been screened in the first infinisee round
'''
from rdkit import Chem
from process_amine_filter_infinisee import get_smiles_from_sdf_or_mol2, extract_molnames_from_sdf


def main():
    sdf_path = 'infinisee/infinisee_patent1/infinisee_patent1_30k.sdf'
    smiles = get_smiles_from_sdf_or_mol2(sdf_path)
    names = extract_molnames_from_sdf(sdf_path)

    already_screened = extract_molnames_from_sdf('infinisee/infinisee_enamine_real_space_similar_to_8gcy.sdf')

    intersection = set(names).intersection(set(already_screened))
    remove_molecules_from_sdf_file(sdf_path, intersection)
    print(len(intersection))

def remove_molecules_from_sdf_file(sdf_path, names_to_remove):
    suppl = Chem.SDMolSupplier(sdf_path)
    mols = [m for m in suppl]
    for m in mols:
        if m.GetProp('_Name') in names_to_remove:
            mols.remove(m)
    w = Chem.SDWriter(sdf_path.replace('.sdf', '_filtered.sdf'))
    for m in mols:
        w.write(m)
    w.close()

if __name__=="__main__":
    main()