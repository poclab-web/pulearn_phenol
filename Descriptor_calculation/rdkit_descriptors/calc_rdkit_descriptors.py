import pandas as pd
from pathlib import Path
from rdkit import Chem
import collections
from rdkit.Chem import Descriptors

def main():
    current_file = Path(__file__).resolve()
    parent_dir = current_file.parent.parent
    ph_smiles_txt_path = parent_dir / 'phenol_smiles_all.txt'

    with open(ph_smiles_txt_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    ph_smiles_lines = [line.strip() for line in lines]

    ph_mols_li = [Chem.MolFromSmiles(smiles) for smiles in ph_smiles_lines]
    MolLogP = [Descriptors.MolLogP(mol) for mol in ph_mols_li]
    MolWt = [Descriptors.MolWt(mol) for mol in ph_mols_li]

    num_N = []
    num_O = []
    num_S = []
    num_F = []
    num_Cl = []
    num_Br = []
    
    for mol in ph_mols_li:
        atoms_list = []
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atoms_list.append(symbol)
            
        atom_dict = collections.Counter(atoms_list)
        
        num_N.append(atom_dict["N"])
        num_O.append(atom_dict["O"])
        num_S.append(atom_dict["S"])
        num_F.append(atom_dict["F"])
        num_Cl.append(atom_dict["Cl"])
        num_Br.append(atom_dict["Br"])
    
    df = pd.DataFrame({'smiles': ph_smiles_lines, 'MolLogP': MolLogP, 'MolWt': MolWt,
                       'num_N': num_N, 'num_O': num_O, 'num_S': num_S,
                       'num_F': num_F, 'num_Cl': num_Cl, 'num_Br': num_Br})
    
    df.to_csv('rdkit_descriptors.csv', index=False)

if __name__ == "__main__":
    main()