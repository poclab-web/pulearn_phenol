import os
import glob
import pandas as pd
from pathlib import Path
import yaml
from rdkit import Chem, rdBase, RDLogger
from rdkit.Chem import AllChem, Descriptors, MolStandardize, PandasTools

RDLogger.DisableLog("rdApp.*")

# Code to output an CSV file summarizing Reagents Group (RG) phenols.

def normalization(df):
    """
    Args:
        df (pandas.DataFrame): DataFrame containing the column "ROMol" with RDKit molecule objects.
    Normalizes the "ROMol" column in the DataFrame by converting each molecule to its InChI representation,
    converting it back to a molecule object, and removing any rows where the InChI conversion fails or results in an empty string.
    Returns:
        pandas.DataFrame: DataFrame with the "ROMol" column containing normalized RDKit molecule objects.
    This function also removes any rows where the InChI conversion fails or results in an empty string.
    """
    inchi = []
    for mol in df["ROMol"]:
        try:
            inchi.append(Chem.MolToInchi(mol))
        except:
            inchi.append("")

    df["InChI"] = inchi
    df.drop(columns="ROMol", inplace=True)
    df.drop(df[df["InChI"] == ""].index, inplace=True)

    df["ROMol"] = df["InChI"].map(Chem.MolFromInchi)
    df.drop(columns="InChI", inplace=True)
    df.dropna(subset=["ROMol"], inplace=True)

    return df

def LargestFragmentChoose(df):
    """
    Chooses the largest fragment from each molecule in the "ROMol" column of the DataFrame.
    If a molecule cannot be processed, it is replaced with an empty string, and those rows are subsequently dropped.
    Args:
        df (pandas.DataFrame): DataFrame containing the column "ROMol" with RDKit molecule objects.
    Normalizes the "ROMol" column by selecting the largest fragment of each molecule.
    Returns:
        pandas.DataFrame: DataFrame with the "ROMol" column containing only the largest fragments of the original molecules.
    This function uses RDKit's LargestFragmentChooser to select the largest fragment from each molecule.
    """
    largest_fragments = []
    for mol in df["ROMol"]:
        try:
            largest_fragments.append(MolStandardize.rdMolStandardize.LargestFragmentChooser().choose(mol))
        except:
            largest_fragments.append("")

    df["ROMol"] = largest_fragments
    df.drop(df[df["ROMol"] == ""].index, inplace=True)

    return df

def exclude_isotopes(df):
    """
    Excludes molecules containing isotopes from the DataFrame.
    Args:
        df (pandas.DataFrame): DataFrame containing the column "ROMol" with RDKit molecule objects.
    This function checks each molecule in the "ROMol" column for the presence of isotopes (e.g., 13C, 15N, 17O, 18O) using a SMARTS pattern.
    Returns:
        pandas.DataFrame: DataFrame with the "ROMol" column containing only molecules that do not contain isotopes.
    This function uses RDKit's substructure matching to identify and exclude isotopes.
    """
    isotopes = Chem.MolFromSmarts("[#1,13c,13C,15n,15N,17o,17O,18o,18O]")

    judge_isotopes = []
    for mol in df["ROMol"]:
        if mol.HasSubstructMatch(isotopes):
            judge_isotopes.append(0)
        else:
            judge_isotopes.append(1)
    df["judge_isotopes"] = judge_isotopes
    df.drop(df[df["judge_isotopes"] == 0].index, inplace=True)

    return df

def Stereolysis(df):
    """
    Converts the "ROMol" column in the DataFrame to SMILES format, removes duplicates, and converts back to RDKit molecule objects.
    This function ensures that the SMILES representation is non-isomeric and removes any duplicate molecules based on their SMILES strings.
    Args:
        df (pandas.DataFrame): DataFrame containing the column "ROMol" with RDKit molecule objects.
    This function first converts each molecule in the "ROMol" column to its SMILES representation, ensuring that the SMILES strings are non-isomeric.
    Returns:
        pandas.DataFrame: DataFrame with the "ROMol" column containing RDKit molecule objects derived from unique SMILES strings.
    This function also removes any duplicate SMILES strings and converts them back to RDKit molecule objects.
    """
    SMILES = []
    for mol in df["ROMol"]:
        SMILES.append(Chem.MolToSmiles(mol, isomericSmiles=False))

    df["SMILES"] = SMILES
    df.drop(columns="ROMol", inplace=True)
    df.drop_duplicates(subset="SMILES", inplace=True)

    df["ROMol"] = df["SMILES"].map(Chem.MolFromSmiles)

    return df

def phenol_scaffold(df):
    """
    Checks if the molecules in the "ROMol" column of the DataFrame contain a phenol substructure.
    This function uses a SMARTS pattern to identify phenol substructures and adds a new column "judge_phenol" to indicate whether each molecule contains the phenol substructure.
    Args:
        df (pandas.DataFrame): DataFrame containing the column "ROMol" with RDKit molecule objects.
    This function checks each molecule in the "ROMol" column for the presence of a phenol substructure using a SMARTS pattern.
    Returns:
        pandas.DataFrame: DataFrame with the "judge_phenol" column indicating whether each molecule contains the phenol substructure.
    This function adds a new column "judge_phenol" to the DataFrame, where 1 indicates the presence of a phenol substructure and 0 indicates its absence.
    """
    phenol_substructure_mol = Chem.MolFromSmarts("[OH1]c1ccccc1")

    judge_phenol = []
    for mol in df["ROMol"]:
        if mol.HasSubstructMatch(phenol_substructure_mol):
            judge_phenol.append(1)
        else:
            judge_phenol.append(0)
    
    df["judge_phenol"] = judge_phenol
    df.drop(df[df["judge_phenol"] == 0].index, inplace=True)

    return df

def count_OH(df):
    """
    Function to count the number of phenolic hydroxyl groups in the phenolic structures.
    This function uses the RDKit Descriptors module to count the number of phenolic hydroxyl groups in each molecule.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'ROMol' that holds RDKit Mol objects.
    Returns:
        pandas.DataFrame: DataFrame with an additional column 'fr_Ar_OH' that contains exactly one phenolic hydroxyl group.
    The 'fr_Ar_OH' column will have a value of 1 if the molecule contains one hydroxyl group, and 0 otherwise.
    """
    df["fr_Ar_OH"] = df["ROMol"].map(Descriptors.fr_Ar_OH)
    df.drop(df[df["fr_Ar_OH"] != 1].index, inplace=True)

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
    phenol_substructure_mol = Chem.MolFromSmarts("[OH1]c1ccccc1")
    replace_core = [AllChem.ReplaceCore(mol, phenol_substructure_mol) for mol in df["ROMol"]]
    replace_core_smiles = [Chem.MolToSmiles(mol) for mol in replace_core]

    judge_naphthol = []
    for smiles1 in replace_core_smiles:
        smiles2 = smiles1.split(".")
        x = len([smiles3 for smiles3 in smiles2 if smiles3.count("*") >= 2])
        if x >= 1:
            judge_naphthol.append(0)
        else:
            judge_naphthol.append(1)

    df["judge_naphthol"] = judge_naphthol
    df.drop(df[df["judge_naphthol"] == 0].index, inplace=True)

    return df

def exclude_ion(df):
    """
    Function to exclude ionized phenolic structures from the DataFrame.
    This function checks if the phenolic structures contain ionized groups and excludes them from the DataFrame.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'SMILES' that holds SMILES representations of phenolic structures.
    Returns:
        pandas.DataFrame: DataFrame with ionized phenolic structures excluded.
    """
    judge_ion = []
    for smiles in df["SMILES"]:
        if "+" in smiles:
            judge_ion.append(0)
        elif "-" in smiles:
            judge_ion.append(0)
        else:
            judge_ion.append(1)
    
    df["judge_ion"] = judge_ion
    df.drop(df[df["judge_ion"] == 0].index, inplace=True)

    return df

def atom_filter(df):
    """
    Function to filter molecules based on the presence of specific atoms.
    This function checks if the molecules in the DataFrame contain only specific target atoms and excludes those that do not.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'ROMol' that holds RDKit Mol objects.
    Returns:
        pandas.DataFrame: DataFrame containing only those molecules that consist solely of the target atoms.
    The target atoms are defined as: ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br'].
    """
    atoms = ["C", "N", "O", "F", "S", "Cl", "Br"]
    
    judge_atoms = []
    for mol in df["ROMol"]:
        li = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in atoms:
                li.append(1)
            else:
                li.append(0)
        
        if sum(li) != mol.GetNumAtoms():
            judge_atoms.append(0)
        else:
            judge_atoms.append(1)
    
    df["judge_atoms"] = judge_atoms
    df.drop(df[df["judge_atoms"] == 0].index, inplace=True)

    return df

def check_ortho(df):
    """
    Function to check if both ortho positions of the phenolic structures are substituted.
    This function uses a SMARTS pattern to check if both ortho positions of the phenolic structures are substituted.
    Args:
        df (pandas.DataFrame): DataFrame containing molecules with a column 'ROMol' that holds RDKit Mol objects.
    Returns:
        pandas.DataFrame: DataFrame with an additional column 'judge_ortho' indicating whether both ortho positions are substituted.
    The 'judge_ortho' column will have a value of 1 if both ortho positions are substituted, and 0 otherwise.
    """
    target = Chem.MolFromSmarts("[OH1]c1c(*)cccc1*")

    judge_ortho = []
    for mol in df["ROMol"]:
        if mol.HasSubstructMatch(target):
            judge_ortho.append(0)
        else:
            judge_ortho.append(1)
    
    df["judge_ortho"] = judge_ortho
    df.drop(df[df["judge_ortho"] == 0].index, inplace=True)

    return df

def main(sdf_folder, result_folder):
    """
    Args:
        sdf_folder (str): Name of the folder containing the SDF file.
        result_folder (str): Name of the folder where the result CSV file will be saved.
    """
    current_file = Path(__file__).resolve()
    parent_dir = current_file.parent.parent
    sdf_path = parent_dir / f'{sdf_folder}'
    result_path = parent_dir / f'{result_folder}'

    sdf_list = glob.glob(f'{sdf_path}/*.sdf')
    df_list=[]
    for sdf in sdf_list:
        df = PandasTools.LoadSDF(sdf)
        df_list.append(df)
    
    df1 = pd.concat(df_list, ignore_index=True)
    df1 = normalization(df1)
    df1 = LargestFragmentChoose(df1)
    df1 = exclude_isotopes(df1)
    df1 = Stereolysis(df1)
    df1 = phenol_scaffold(df1)
    df1 = count_OH(df1)
    df1 = exclude_naphthol(df1)
    df1 = exclude_ion(df1)
    df1 = atom_filter(df1)
    df1 = check_ortho(df1)

    df1["inchi"] = df1["ROMol"].map(Chem.MolToInchi)
    df1["inchikey"] = df1["ROMol"].map(Chem.MolToInchiKey)
    df1['smiles'] = df1['ROMol'].map(Chem.MolToSmiles)

    df1 = df1[['inchi', 'inchikey', 'smiles']]
    df1.dropna(inplace=True)
    df1.drop_duplicates(subset="inchi", inplace=True)

    df1.to_csv(f'{result_path}/phenols_RG.csv', index=False)

if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))
    with open("Preprocessing_RG.yaml", mode="rb") as yml:
        settings = yaml.safe_load(yml)
    
    sdf_folder = settings["sdf_folder"]
    result_folder = settings["result_folder"]

    main(sdf_folder, result_folder)