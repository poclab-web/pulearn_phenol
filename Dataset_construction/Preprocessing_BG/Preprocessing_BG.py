import os
import sys
import pandas as pd
from pathlib import Path
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, PandasTools, Descriptors

# Code to output an Excel file summarizing Biological Group (BG) phenols.

def biphenol_scaffold(df):
    """
    Function to filter the DataFrame for biphenol structures.
    This function checks if the molecules in the DataFrame contain a biphenol scaffold.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'mol_pro' that holds RDKit Mol objects.
    Returns:
        pandas.DataFrame: Filtered DataFrame containing only those molecules that match the biphenol scaffold.
    The biphenol scaffold is defined by the SMARTS pattern '[OH]c1c(-c2c([OH])cccc2)cccc1'.
    """

    biphenol_ = '[OH]c1c(-c2c([OH])cccc2)cccc1'
    biphenol_mol = Chem.MolFromSmarts(biphenol_)

    judge = []
    for mol in df['mol_pro']:
        if mol.HasSubstructMatch(biphenol_mol):
            judge.append(1)
        else:
            judge.append(0)
    df['judge_biphenol'] = judge
    df = df[df['judge_biphenol']==1]
    
    return df

def Disassembly(df):
    """
    Function to disassemble biphenol structures into two phenolic structures.
    This function uses a SMARTS reaction to break down biphenol structures into their constituent parts.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'mol_pro' that holds RDKit Mol objects.
    The function applies a reaction to each molecule in 'mol_pro' to generate two phenolic products.
    Returns:
        pandas.DataFrame: DataFrame with two new columns 'phenol_smiles' and 'phenol_smiles2',
        which contain the SMILES representations of the two phenolic products obtained from the disassembly.
    The reaction is defined by the SMARTS pattern:
    '[O:1]-[c:2]1:[c:3](-[c:10]2:[c:11]:[c:12]:[c:13]:[c:14]:[c:9]:2-[O:8]):[c:4]:[c:5]:[c:6]:[c:7]:1 >> [O:1][c:2]1[c:3][c:4][c:5][c:6][c:7]1.[O:8][c:9]1[c:10][c:11][c:12][c:13][c:14]1'.
    """
    
    smarts = '[O:1]-[c:2]1:[c:3](-[c:10]2:[c:11]:[c:12]:[c:13]:[c:14]:[c:9]:2-[O:8]):[c:4]:[c:5]:[c:6]:[c:7]:1 >> [O:1][c:2]1[c:3][c:4][c:5][c:6][c:7]1.[O:8][c:9]1[c:10][c:11][c:12][c:13][c:14]1'
    reaction = AllChem.ReactionFromSmarts(smarts)
    phenol_smiles = []
    phenol_smiles2 = []
    
    for mol in df['mol_pro']:
        products = reaction.RunReactants([mol])
        try:
            phenol_smiles.append(Chem.MolToSmiles(products[0][0], isomericSmiles=False))
        except:
            phenol_smiles.append(None)
        try:
            phenol_smiles2.append(Chem.MolToSmiles(products[0][1], isomericSmiles=False))
        except:
            phenol_smiles2.append(None)
        
    df['phenol_smiles'] = phenol_smiles
    df['phenol_smiles2'] = phenol_smiles2
        
    df = df.dropna(subset=['phenol_smiles', 'phenol_smiles2'])
    df = df[df['phenol_smiles']==df['phenol_smiles2']]
    
    return df

def count_OH(df):
    """
    Function to count the number of phenolic hydroxyl groups in the phenolic structures.
    This function uses the RDKit Descriptors module to count the number of phenolic hydroxyl groups in each molecule.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'mol_rea' that holds RDKit Mol objects.
    Returns:
        pandas.DataFrame: DataFrame with an additional column 'fr_Ar_OH' that contains exactly one phenolic hydroxyl group.
    The 'fr_Ar_OH' column will have a value of 1 if the molecule contains one hydroxyl group, and 0 otherwise.
    """
    
    df['fr_Ar_OH'] = df["mol_rea"].map(Descriptors.fr_Ar_OH)
    df = df[df['fr_Ar_OH']==1]
    
    return df

def exclude_naphthol(df):
    """
    Function to exclude condensed polycyclic structures such as naphthols from the DataFrame.
    This function checks if the phenolic structures contain condensed polycyclic structures and excludes them from the DataFrame.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'mol_rea' that holds RDKit Mol objects.
    Returns:
        pandas.DataFrame: DataFrame with condensed polycyclic structures excluded.
    """
    
    phenol_smarts = '[OH1]c1ccccc1'
    phenol_mol = Chem.MolFromSmarts(phenol_smarts)
    
    replace_core = [AllChem.ReplaceCore(mol, phenol_mol) for mol in df['mol_rea']]
    replace_core_smiles = [Chem.MolToSmiles(mol) for mol in replace_core]
    
    judge_naphthol=[]
    
    for smiles in replace_core_smiles:
        
        smiles2 = smiles.split(".")
        x = len([smiles for smiles in smiles2 if smiles.count('*')>=2])
        
        if x>=1:
            judge_naphthol.append(1)
        else:
            judge_naphthol.append(0)
    
    df['judge_naphthol'] = judge_naphthol
    df = df[df['judge_naphthol']==0]
    
    return df

def exclude_ion(df):
    """
    Function to exclude ionized phenolic structures from the DataFrame.
    This function checks if the phenolic structures contain ionized groups and excludes them from the DataFrame.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'phenol_smiles' that holds SMILES representations of phenolic structures.
    Returns:
        pandas.DataFrame: DataFrame with ionized phenolic structures excluded.
    """
    
    judge_ion=[]
    
    for smiles in df['phenol_smiles']:
        
        if '+' in smiles:
            judge_ion.append(1)
        elif '-' in smiles:
            judge_ion.append(1)
        else:
            judge_ion.append(0)
    
    df['judge_ion'] = judge_ion
    df = df[df['judge_ion']==0]
    
    return df

def MolWt_check(df):
    """
    This function calculates the molecular weights of phenol and biphenol structures, and checks whether the molecular weight of the biphenol is approximately twice that of the phenol.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with columns 'mol_rea' and 'mol_pro' that hold RDKit Mol objects.
    Returns:
        pandas.DataFrame: DataFrame containing molecules that meet the conditions described above.
    """

    df["MolWt_rea"] = df["mol_rea"].map(Descriptors.MolWt)
    df["MolWt_pro"] = df["mol_pro"].map(Descriptors.MolWt)
    
    df['MolWt_check'] = df['MolWt_pro']-df['MolWt_rea']*2+2
    df = df[(df['MolWt_check']>-1)&(df['MolWt_check']<1)]
    
    return df

def atom_filter(df):
    """
    Function to filter molecules based on the presence of specific atoms.
    This function checks if the molecules in the DataFrame contain only specific target atoms and excludes those that do not.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'mol_rea' that holds RDKit Mol objects.
    Returns:
        pandas.DataFrame: DataFrame containing only those molecules that consist solely of the target atoms.
    The target atoms are defined as: ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br'].
    """
    
    target_atoms = ['C','N','O','S','F','Cl','Br']
    delete = []
    
    for mol in df['mol_rea']:
        
        include = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in target_atoms:
                include.append(1)
            else:
                include.append(0)
        
        if sum(include) != mol.GetNumAtoms():
            delete.append(1)
        else:
            delete.append(0)
    
    df['delete'] = delete
    df = df[df['delete']==0]
    
    return df

def Stereolysis(df):
    """
    Function to remove duplicate molecules in a DataFrame.
    This function checks the duplicate molecules in a DataFrame by InChI.
    Args: 
        df (pandas.DataFrame): DataFrame containing molecules with a column 'mol_rea' that holds RDKit Mol objects.
    Returns:
        pandas.DataFrame: DataFrame containing only unique molecules based on their InChI representation.
    """
    
    df_2 = df.copy()
    df_2['inchi_rea'] = [Chem.MolToInchi(mol) for mol in df_2['mol_rea']]
    df_2['inchikey_rea'] = [Chem.MolToInchiKey(mol) for mol in df_2['mol_rea']]
    df_2 = df_2[['inchi_rea', 'inchikey_rea', 'phenol_smiles', 'smiles_pro']]
    df_2 = df_2.rename(columns={'phenol_smiles':'smiles_rea'})
    df_2 = df_2.drop_duplicates(subset='inchi_rea')
    
    return df_2

def main(sdf_folder, result_folder):
    """
    Args:
        sdf_folder (str): Name of the folder containing the SDF file.
        result_folder (str): Name of the folder where the result Excel file will be saved.
    """

    current_file = Path(__file__).resolve()
    parent_dir = current_file.parent.parent
    sdf_path = parent_dir / f'{sdf_folder}'
    result_path = parent_dir / f'{result_folder}'

    sdf_list = list(sdf_path.rglob('*.sdf'))
    df_list=[]
    for sdf in sdf_list:
        df = PandasTools.LoadSDF(sdf)
        df_list.append(df)
    
    df = pd.concat(df_list, ignore_index=True)
    df["smiles_pro"] = [Chem.MolToSmiles(mol, isomericSmiles=False) for mol in df['ROMol']]
    df['mol_pro'] = [Chem.MolFromSmiles(smiles) for smiles in df['smiles_pro']]
    df = biphenol_scaffold(df)
    df = Disassembly(df)
    df['mol_rea'] = [Chem.MolFromSmiles(smiles) for smiles in df['phenol_smiles']]
    df = count_OH(df)
    df = exclude_naphthol(df)
    df = exclude_ion(df)
    df = MolWt_check(df)
    df = atom_filter(df)
    df = Stereolysis(df)
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol="smiles_rea", molCol='image_rea')
    PandasTools.SaveXlsxFromFrame(df, f'{result_path}/phenols_BG.xlsx',
                                  molCol='image_rea', size = (200,200))

if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))
    with open("Preprocessing_BG.yaml", mode="rb") as yml:
        settings = yaml.safe_load(yml)
    
    sdf_folder = settings["sdf_folder"]
    result_folder = settings["result_folder"]

    main(sdf_folder, result_folder)