import os
import pandas as pd
import yaml
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdmolops, AllChem, Draw

def get_indexes_SG_BG(smiles_ph, smiles_bi):
    """Get the index of atoms that compose the phenol skeleton.

    - Returns a list of RDKit atom indexes in the order of
      O, next_O, o_1, m_1, p, m_2, o_2.
    - For the Synthesized Group, o_1 was defined as the ortho position of the monomeric phenol 
      that formed a bond with its counterpart in oxidative homocoupling dimers retrieved from SciFinder.
    - For the Biological Group, o_1 was similarly defined as the ortho position of the monomeric phenol
      that bonded with its counterpart in dimers obtained from databases such as SciFinder, ChEMBL, ChEBI, and PubChem.
    - In the case of different substituents on both sides of the meta position, the index was entered 
      by looking at the images output in the folder 'images_ck'.

    Args:
        smiles_ph (str): SMILES of phenols.
        smiles_bi (str): SMILES of biphenols.

    Returns:
        list: [O, next_O, o_1, m_1, p, m_2, o_2]
    """

    ph_substructure = '[OH1]c1ccccc1'
    ph_substructure_mol = Chem.MolFromSmarts(ph_substructure)

    mol = Chem.MolFromSmiles(smiles_ph)
    match_atoms = list(mol.GetSubstructMatches(ph_substructure_mol))
    match_atoms = list(match_atoms[0])

    for atom in mol.GetAtoms():
        if (atom.GetIdx() in match_atoms) & (atom.GetSymbol() == 'O'):
            O = atom.GetIdx()
            match_atoms.remove(atom.GetIdx())
            break
    
    subsituents = isolating_substituents_SG_BG(smiles_bi)

    if (subsituents[0] == subsituents[2]) | (subsituents[3] != '[2H][H]'):
        for atom in mol.GetAtoms():
            if atom.GetIdx() != O:
                dis_from_O = rdmolops.GetShortestPath(mol, O, atom.GetIdx())

                if (len(dis_from_O) == 3) & (len(atom.GetNeighbors()) == 2):
                    o_1 = atom.GetIdx()
                elif len(dis_from_O) == 2:
                    next_O = atom.GetIdx()          
        match_atoms.remove(next_O)
        match_atoms.remove(o_1)

        dis_from_O = rdmolops.GetShortestPath(mol, O, o_1)

        for ph_atom in match_atoms:
            dis_from_start_C = rdmolops.GetShortestPath(mol, o_1, ph_atom)
            and_list = set(dis_from_O) & set(dis_from_start_C)

            if (len(dis_from_start_C) == 2) & (len(and_list) == 1):
                m_1 = ph_atom
            elif (len(dis_from_start_C) == 3) & (len(and_list) == 1):
                p = ph_atom
            elif len(dis_from_start_C) == 4:
                m_2 = ph_atom
            elif (len(dis_from_start_C) == 3) & (len(and_list) == 2):
                o_2 = ph_atom
    
    elif not (subsituents[0] != '[2H][H]') & (subsituents[2] != '[2H][H]'):
        for atom in mol.GetAtoms():
            if atom.GetIdx() != O:
                dis_from_O = rdmolops.GetShortestPath(mol, O, atom.GetIdx())
                
                if (len(dis_from_O) == 4) & (len(atom.GetNeighbors()) == 2):
                        start_C = atom.GetIdx()
                
                elif len(dis_from_O) == 2:
                    next_O = atom.GetIdx()

        match_atoms.remove(start_C)
        match_atoms.remove(next_O)
        dis_from_O = rdmolops.GetShortestPath(mol, O, start_C)

        if subsituents[0] == '[2H][H]':
            m_1 = start_C

            for ph_atom in match_atoms:
                dis_from_start_C = rdmolops.GetShortestPath(mol, start_C, ph_atom)
                and_list = set(dis_from_O) & set(dis_from_start_C)

                if (len(dis_from_start_C) == 2) & (len(and_list) == 2):
                    o_1 = ph_atom
                elif (len(dis_from_start_C) == 2) & (len(and_list) == 1):
                    p = ph_atom
                elif (len(dis_from_start_C) == 3) & (len(and_list) == 1):
                    m_2 = ph_atom
                elif (len(dis_from_start_C) == 4):
                    o_2 = ph_atom
        
        elif subsituents[2] == '[2H][H]':
            m_2 = start_C

            for ph_atom in match_atoms:
                dis_from_start_C = rdmolops.GetShortestPath(mol, start_C, ph_atom)
                and_list = set(dis_from_O) & set(dis_from_start_C)

                if (len(dis_from_start_C) == 2) & (len(and_list) == 2):
                    o_2 = ph_atom
                elif (len(dis_from_start_C) == 2) & (len(and_list) == 1):
                    p = ph_atom
                elif (len(dis_from_start_C) == 3) & (len(and_list) == 1):
                    m_1 = ph_atom
                elif len(dis_from_start_C) == 4:
                    o_1 = ph_atom
    
    else:
        for atom in mol.GetAtoms():
            if atom.GetIdx() != O:
                dis_from_O = rdmolops.GetShortestPath(mol, O, atom.GetIdx())

                if len(dis_from_O)==2:
                    next_O = atom.GetIdx()
            
            atom.SetProp("atomLabel", str(atom.GetIdx()))
        
        mol_bi = Chem.MolFromSmiles(smiles_bi)
        
        os.makedirs('images_ck', exist_ok=True)
        Draw.MolToFile(mol, f'images_ck/{smiles_ph}.png')
        Draw.MolToFile(mol_bi, f'images_ck/{smiles_bi}.png')

        print(f'Please fill in the following fields with reference to images_ck/{smiles_ph}.png and images_ck/{smiles_bi}.png')
        o_1 = int(input("Please enter the atom index for o_1 :")) 
        m_1 = int(input("Please enter the atom index for m_1 :")) 
        p = int(input("Please enter the atom index for p :")) 
        m_2 = int(input("Please enter the atom index for m_2 :"))
        o_2 = int(input("Please enter the atom index for o_2 :"))
    
    return [O, next_O, o_1, m_1, p, m_2, o_2]

def get_indexes_RG(smiles):
    """Get the index of atoms that compose the phenol skeleton.

    - Returns a list of RDKit atom indexes in the order of
      O, next_O, o_1, m_1, p, m_2, o_2.
    - The higher Cahn-Ingold-Prelog (CIP) priority at the ortho position
      is designated as o_2.
    - If the CIP priority at the ortho positions are equal,
      the higher CIP priority at the meta position is designated as m_2.

    Args:
        smiles (str): SMILES of phenols reagents.

    Returns:
        list: [O, next_O, o_1, m_1, p, m_2, o_2]
    """
    mol = Chem.MolFromSmiles(smiles)
    phenol_substructure = Chem.MolFromSmarts("[OH1]c1ccccc1")
    match_atoms = mol.GetSubstructMatch(phenol_substructure)

    C = []
    for atom in mol.GetAtoms():
        if atom.GetIdx() in match_atoms:
            if atom.GetSymbol() == "O":
                O = atom.GetIdx()
            else:
                C.append(atom.GetIdx())

    ortho = []
    meta = []
    for i in C:
        distance_from_O = rdmolops.GetShortestPath(mol, O, i)
        if len(distance_from_O) == 2:
            next_O = i
        elif len(distance_from_O) == 3:
            ortho.append(i)
        elif len(distance_from_O) == 4:
            meta.append(i)
        elif len(distance_from_O) == 5:
            p = i

    ortho_cip = []
    for i in ortho:
        ortho_cip.append(int(mol.GetAtomWithIdx(i).GetProp("_CIPRank")))

    if ortho_cip[0] < ortho_cip[1]:
        o_1 = ortho[0]
        o_2 = ortho[1]
        for i in meta:
            distance_from_o_1 = rdmolops.GetShortestPath(mol, o_1, i)
            if len(distance_from_o_1) == 2:
                m_1 = i
            elif len(distance_from_o_1) == 4:
                m_2 = i

    elif ortho_cip[0] > ortho_cip[1]:
        o_1 = ortho[1]
        o_2 = ortho[0]
        for i in meta:
            distance_from_o_1 = rdmolops.GetShortestPath(mol, o_1, i)
            if len(distance_from_o_1) == 2:
                m_1 = i
            elif len(distance_from_o_1) == 4:
                m_2 = i

    else:
        meta_cip = []
        for i in meta:
            meta_cip.append(int(mol.GetAtomWithIdx(i).GetProp("_CIPRank")))

        if meta_cip[0] < meta_cip[1]:
            m_1 = meta[0]
            m_2 = meta[1]
            for i in ortho:
                distance_from_m_1 = rdmolops.GetShortestPath(mol, m_1, i)
                if len(distance_from_m_1) == 2:
                    o_1 = i
                elif len(distance_from_m_1) == 4:
                    o_2 = i

        elif meta_cip[0] > meta_cip[1]:
            m_1 = meta[1]
            m_2 = meta[0]
            for i in ortho:
                distance_from_m_1 = rdmolops.GetShortestPath(mol, m_1, i)
                if len(distance_from_m_1) == 2:
                    o_1 = i
                elif len(distance_from_m_1) == 4:
                    o_2 = i

        else:
            o_1 = ortho[0]
            o_2 = ortho[1]
            for i in meta:
                distance_from_o_1 = rdmolops.GetShortestPath(mol, o_1, i)
                if len(distance_from_o_1) == 2:
                    m_1 = i
                elif len(distance_from_o_1) == 4:
                    m_2 = i

    return [O, next_O, o_1, m_1, p, m_2, o_2]

def isolating_substituents_SG_BG(smiles):
    """Isolate substituents.

    - Returns a list of SMILES of substituents in the order of m_1, p, m_2, o.
    - 2H (deuterium) at the end of the substituents.

    Args:
        smiles (str): SMILES of biphenols.

    Returns:
        list: [m_1, p, m_2, o]
    """
    
    ph_substructure = '[OH1]c1ccccc1'
    ph_substructure_mol = Chem.MolFromSmarts(ph_substructure)

    mol_ph = Chem.MolFromSmiles(smiles)

    replace_side = AllChem.ReplaceSidechains(mol_ph, ph_substructure_mol)
    replace_core = AllChem.ReplaceCore(mol_ph, ph_substructure_mol)

    replace_side_smiles = Chem.MolToSmiles(replace_side)
    replace_core_smiles = Chem.MolToSmiles(replace_core)

    metal = ['[Li]','[Na]','[K]','[Rb]','[Cs]']

    for i in range(5):
        if f'[{i+1}*]' in replace_side_smiles:
            replace_side_smiles = replace_side_smiles.replace(f'[{i+1}*]', metal[i])
        
        if f'[{i+1}*]' in replace_core_smiles:
            replace_core_smiles = replace_core_smiles.replace(f'[{i+1}*]', metal[i])
    
    replace_core_smiles_li = replace_core_smiles.split('.')
    replace_side_2 = Chem.MolFromSmiles(replace_side_smiles)

    replace_core_onePhenol = [Chem.MolFromSmiles(smiles) for smiles in replace_core_smiles_li]

    for sub in replace_core_onePhenol:
        if sub.HasSubstructMatch(ph_substructure_mol):
            ph_sub = replace_core_smiles_li[replace_core_onePhenol.index(sub)]
    
    for metal_atom in metal:
        if metal_atom in ph_sub:
            for atom in replace_side_2.GetAtoms():
                if atom.GetSymbol() == metal_atom.strip("[]"):
                    start_atom_Idx = atom.GetIdx()
    
    # Get phenol oxygen index (a) and metal atom (stored in b)
    b = []
    for atom in replace_side_2.GetAtoms():
        if atom.GetSymbol() == 'O':
            a = atom.GetIdx()
        
        elif (atom.GetSymbol() != 'C') & (atom.GetIdx() != start_atom_Idx):
            b.append(atom)
    
    dis_from_O = rdmolops.GetShortestPath(replace_side_2, a, start_atom_Idx)

    m_1 = []
    p = []
    m_2 = []
    o = []

    for i in range(len(b)):
        dis_from_reaction_point = rdmolops.GetShortestPath(replace_side_2, start_atom_Idx, b[i].GetIdx())
        and_list = set(dis_from_O) & set(dis_from_reaction_point)
        
        smiles = list(filter(lambda x: f'{b[i].GetSymbol()}' in x, replace_core_smiles_li))[0]
        
        if (len(dis_from_reaction_point) == 4) & (len(and_list) == 2):
            m_1.append(smiles.replace(f'[{b[i].GetSymbol()}]','[2H]'))
        elif (len(dis_from_reaction_point) == 5) & (len(and_list) == 2):
            p.append(smiles.replace(f'[{b[i].GetSymbol()}]','[2H]'))
        elif (len(dis_from_reaction_point) == 6):
            m_2.append(smiles.replace(f'[{b[i].GetSymbol()}]','[2H]'))
        elif (len(dis_from_reaction_point) == 5) & (len(and_list) == 3):
            o.append(smiles.replace(f'[{b[i].GetSymbol()}]','[2H]'))
        
    if not m_1:
        m_1.append('[2H][H]')
    if not p:
        p.append('[2H][H]')
    if not m_2:
        m_2.append('[2H][H]')
    if not o:
        o.append('[2H][H]')
    
    li = [m_1[0], p[0], m_2[0], o[0]]

    return li

def isolating_substituents_RG(smiles):
    """Isolate substituents.

    - Returns a list of SMILES of substituents in the order of m_1, p, m_2, o.
    - 2H (deuterium) at the end of the substituents.

    Args:
        smiles (str): SMILES of phenols reagents.

    Returns:
        list: [m_1, p, m_2, o]
    """
    mol = Chem.MolFromSmiles(smiles)
    indexes = get_indexes_RG(smiles)

    li = []
    for i in range(3, 7):
        substituent = []

        fragments = Chem.FragmentOnBonds(mol, [indexes[i]])
        fragments_smiles = Chem.MolToSmiles(fragments)
        list_fragments_smiles = fragments_smiles.split(".")

        if len(list_fragments_smiles) == 2:
            for smiles in list_fragments_smiles:
                if f"[{indexes[i]}*]" in smiles:
                    substituent.append(smiles.replace(f"{indexes[i]}*", "2H"))
                    break

        if len(substituent) == 0:
            fragments = Chem.FragmentOnBonds(mol, [indexes[i]-1])
            fragments_smiles = Chem.MolToSmiles(fragments)
            list_fragments_smiles = fragments_smiles.split(".")

            if len(list_fragments_smiles) == 2:
                for smiles in list_fragments_smiles:
                    if f"[{indexes[i]}*]" in smiles:
                        substituent.append(smiles.replace(f"{indexes[i]}*", "2H"))
                        break

        if len(substituent) == 0:
            substituent.append("[2H][H]")

        li.append(substituent[0])

    return li

def main(ph_smiles_txt, bi_smiles_txt, keyword):
    current_file = Path(__file__).resolve()
    parent_dir = current_file.parent.parent
    ph_smiles_txt_path = parent_dir / f'{ph_smiles_txt}'
    bi_smiles_txt_path = parent_dir / f'{bi_smiles_txt}'

    with open(ph_smiles_txt_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    ph_smiles_lines = [line.strip() for line in lines]
    
    try:
        with open(bi_smiles_txt_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        bi_smiles_lines = [line.strip() for line in lines]
    except:
        pass
    
    subs = []
    indexes = []
    if keyword == 'SG_BG':
        for ph_smiles, bi_smiles in zip(ph_smiles_lines, bi_smiles_lines):
            subs.append(isolating_substituents_SG_BG(bi_smiles))
            indexes.append(get_indexes_SG_BG(ph_smiles, bi_smiles))
        
        df_subs = pd.DataFrame(subs, columns=['m_1', 'p', 'm_2', 'o'])
        df_indexes = pd.DataFrame(indexes, columns=['O', 'next_O', 'o_1', 'm_1', 'p', 'm_2', 'o_2'])
        df_subs.insert(0, 'smiles', ph_smiles_lines)
        df_indexes.insert(0, 'smiles', ph_smiles_lines)

        df_subs.to_csv(f'isolating_substituents_{keyword}.csv', index=False)
        df_indexes.to_csv(f'get_indexes_{keyword}.csv', index=False)
    
    elif keyword == 'RG':
        for ph_smiles in ph_smiles_lines:
            subs.append(isolating_substituents_RG(ph_smiles))
            indexes.append(get_indexes_RG(ph_smiles))
        
        df_subs = pd.DataFrame(subs, columns=['m_1', 'p', 'm_2', 'o'])
        df_indexes = pd.DataFrame(indexes, columns=['O', 'next_O', 'o_1', 'm_1', 'p', 'm_2', 'o_2'])
        df_subs.insert(0, 'smiles', ph_smiles_lines)
        df_indexes.insert(0, 'smiles', ph_smiles_lines)

        df_subs.to_csv(f'isolating_substituents_{keyword}.csv', index=False)
        df_indexes.to_csv(f'get_indexes_{keyword}.csv', index=False)
    
    else:
        print('The keyword you entered is not correct')

if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))
    with open("tools.yaml", mode="rb") as yml:
        settings = yaml.safe_load(yml)
    
    ph_smiles_txt = settings["ph_smiles_txt"]
    bi_smiles_txt = settings["bi_smiles_txt"]
    keyword = settings["keyword"]

    main(ph_smiles_txt, bi_smiles_txt, keyword)