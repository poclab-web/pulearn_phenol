import os
import pandas as pd
from pathlib import Path
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, PandasTools, Descriptors

# Code to output an Excel file summarizing Synthesized Group (SG) phenols.

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
    ph_substructure = '[OH1]c1ccccc1'
    ph_substructure_mol = Chem.MolFromSmarts(ph_substructure)
    
    replace_core = [AllChem.ReplaceCore(mol,ph_substructure_mol) for mol in df['mol_rea']]
    replace_core_smiles = [Chem.MolToSmiles(mol) for mol in replace_core]
    
    judge_naphthol = []
    
    for smiles in replace_core_smiles:
        
        smiles2 = smiles.split(".")
        x = len([smiles for smiles in smiles2 if smiles.count('*')>=2])
        
        if x>=1:
            judge_naphthol.append(0)
        else:
            judge_naphthol.append(1)
    
    df['judge_naphthol'] = judge_naphthol
    df = df[df['judge_naphthol']==1]
    
    return df

def exclude_ion(df):
    """
    Function to exclude ionized phenolic structures from the DataFrame.
    This function checks if the phenolic structures contain ionized groups and excludes them from the DataFrame.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'smiles_rea' that holds SMILES representations of phenolic structures.
    Returns:
        pandas.DataFrame: DataFrame with ionized phenolic structures excluded.
    """
    judge_ion = []
    
    for smiles in df['smiles_rea']:
        if '+' in smiles:
            judge_ion.append(0)
        elif '-' in smiles:
            judge_ion.append(0)
        else:
            judge_ion.append(1)
    
    df['judge_ion'] = judge_ion
    df = df[df['judge_ion']==1]
    
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
    atoms = ['C','N','O','S','F','Cl','Br']
    delete = []
    
    for mol in df['mol_rea']:
        include = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in atoms:
                include.append(1)
            else:
                include.append(0)
        
        if sum(include) != mol.GetNumAtoms():
            delete.append(1)
        else:
            delete.append(0)
    
    df['delete'] = delete
    df = df[df['delete']==0]
    df = df.drop('delete',axis=1)
    
    return df

def Stereolysis(df):
    """
    Function to remove duplicate molecules in a DataFrame.
    This function checks the duplicate molecules in a DataFrame by InChI.
    Args: 
        df (pandas.DataFrame): DataFrame containing molecules with a column 'mol_rea' and 'mol_pro' that holds RDKit Mol objects.
    Returns:
        pandas.DataFrame: DataFrame containing only unique molecules based on their InChI representation.
    """
    df['smiles_rea'] = [Chem.MolToSmiles(mol, isomericSmiles=False) for mol in df['mol_rea']]
    df['smiles_pro'] = [Chem.MolToSmiles(mol, isomericSmiles=False) for mol in df['mol_pro']]
    df['mol_rea'] = [Chem.MolFromSmiles(smiles) for smiles in df['smiles_rea']]
    df['inchi_rea'] = [Chem.MolToInchi(mol) for mol in df['mol_rea']]
    df['inchikey_rea'] = [Chem.MolToInchiKey(mol) for mol in df['mol_rea']]
    
    df = df[['inchi_rea','inchikey_rea','smiles_rea','smiles_pro', 'mol_rea',
             'RGT','SOL','CAT','ylds','reference']]
    
    return df

def main(after_rdf_analyze_file, result_folder):
    """
    Args:
        after_rdf_analyze_file (str): Name of the Excel file output by rdf_analyze.py.
        result_folder (str): Name of the folder where the result Excel file will be saved.
    """
    current_file = Path(__file__).resolve()
    parent_dir = current_file.parent.parent
    file_path = parent_dir / f'output/{after_rdf_analyze_file}'
    result_path = parent_dir / f'{result_folder}'

    df = pd.read_excel(f'{file_path}')
    df['mol_rea'] = [Chem.MolFromSmiles(smiles) for smiles in df['smiles_rea']]
    df['mol_pro'] = [Chem.MolFromSmiles(smiles) for smiles in df['smiles_pro']]
    df = df.query('num_rea == 1')
    df = df.query('mol_rea.notnull()', engine='python')
    df = df.query('mol_pro.notnull()', engine='python')

    df = MolWt_check(df)
    df = count_OH(df)
    df = exclude_naphthol(df)
    df = exclude_ion(df)
    df = atom_filter(df)
    df = Stereolysis(df)

    df2 = df.loc[:, :'mol_rea']
    df2 = df2.drop_duplicates(subset='inchi_rea')
    PandasTools.SaveXlsxFromFrame(df2, f'{result_path}/phenols_SG.xlsx', 
                                  molCol='mol_rea', size=(200, 200))
    
    df = df.drop('mol_rea', axis=1)
    df.to_csv(f'{result_path}/reactions_SG.csv', index=False)

if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))
    with open("Preprocessing_SG.yaml", mode="rb") as yml:
        settings = yaml.safe_load(yml)
    
    after_rdf_analyze_file = settings["after_rdf_analyze_file"]
    result_folder = settings["result_folder"]

    main(after_rdf_analyze_file, result_folder)