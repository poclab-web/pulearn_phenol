import os
import glob
import pandas as pd
import numpy as np
from pathlib import Path
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, PandasTools, Descriptors

def read_rdf(rdf_list):
    """
    Read RDF files and extract reaction strings.
    This function reads multiple RDF files, extracts the reaction strings,
    and returns a DataFrame containing these strings along with their corresponding reactions.
    Args:
        rdf_list (list): List of RDF file paths to be processed.
    The RDF files are expected to be in a specific format where reaction strings can be identified.
    Returns:
        pandas.DataFrame: A DataFrame containing the reaction strings and their corresponding RDKit reaction objects.
    """
    list_ = []
    for _ in rdf_list:
        with open(_, encoding="utf_8") as f:
            rdf = f.read()
        rdf_n = rdf.split("\n")
        rdf_RXN = '\n'.join(rdf_n[3:]).split("$RXN")
        rxn_str = ["$RXN"+s for s in rdf_RXN[1:]]
        df = pd.DataFrame(rxn_str,columns=['rxn_str'])
        list_.append(df)
    df = pd.concat(list_, ignore_index=True, axis=0)
    
    rxns=[]
    for i in range(len(df)):
        with open('rxn.txt', 'w', encoding='utf_8') as f:
            f.write(df["rxn_str"][i])
        rxn = AllChem.ReactionFromRxnFile("rxn.txt")
        rxns.append(rxn)
    df["rxn"] = rxns
    return df

def rxn_analyze(df):
    """
    Analyze the reactions in the DataFrame to extract information about reactants and products.
    This function adds columns to the DataFrame for the number of reactants and products,
    as well as the reactants and products themselves.
    Args:
        df (pandas.DataFrame): DataFrame containing reaction strings and RDKit reaction objects.
    The DataFrame is expected to have a column 'rxn' containing RDKit reaction objects.
    The function will add the following columns to the DataFrame:
        - 'num_rea': Number of reactants in each reaction.
        - 'num_pro': Number of products in each reaction.
        - 'mols_rea': List of reactant molecules for each reaction.
        - 'mols_pro': List of product molecules for each reaction.
    Returns:
        pandas.DataFrame: DataFrame with additional columns for reactants and products.
    """
    df['num_rea'] = [rxn.GetNumReactantTemplates() for rxn in df['rxn']]
    df['num_pro'] = [rxn.GetNumProductTemplates() for rxn in df['rxn']]
    df['mols_rea'] = [rxn.GetReactants() for rxn in df['rxn']]
    df['mols_pro'] = [rxn.GetProducts() for rxn in df['rxn']]
    return df

def get_conditions(df):
    """
    Extracts reaction conditions from the DataFrame and adds them as new columns.
    This function processes the reaction strings in the DataFrame to extract reaction conditions such as reagents, solvent, and catalyst.
    Args:
        df (pandas.DataFrame): DataFrame containing reaction strings.
    The DataFrame is expected to have a column 'rxn_str' containing reaction strings in a specific format.
    Returns:
        pandas.DataFrame: DataFrame with additional columns for reagents, solvent, and catalyst.
    """
    for name_ in ['RGT','SOL','CAT']:
        pro = []
        for i in range(len(df)):
            rxn_n = df['rxn_str'][i].split('\n')
            pro_ = []
            i=1
            while True:
                try:
                    pro_.append(rxn_n[rxn_n.index('$DTYPE RXN:VAR(1):'+name_+'('+str(i)+'):CAS_RN')+1].strip('$DATUM '))
                    i+=1
                except: break
            pro.append(pro_)
        df[name_] = pro
    return df

def get_smiles_rea(df):
    """
    Converts the reactant molecules in the DataFrame to SMILES format.
    This function processes the reactant molecules in the DataFrame and converts them to SMILES strings.
    Args:
        df (pandas.DataFrame): DataFrame containing reaction strings and RDKit reaction objects.
    The DataFrame is expected to have a column 'mols_rea' containing lists of RDKit molecule objects for reactants.
    Returns:
        pandas.DataFrame: DataFrame with an additional column 'smiles_rea' containing SMILES strings of the reactants.
    The function handles cases where the conversion to SMILES fails by appending 'error' to the list.
    """
    smiles_reas = []
    for i in range(len(df)):
        try:
            mol_rea = df['mols_rea'][i][0]
            smiles_rea = Chem.MolToSmiles(mol_rea)
            smiles_reas.append(smiles_rea)
        except:
            smiles_reas.append('error')
    df['smiles_rea'] = smiles_reas
    return df

def get_smiles_pro_yld(df):
    """
    Extracts the product SMILES and yields from the DataFrame.
    This function processes the product molecules in the DataFrame and extracts their SMILES representations.
    It also retrieves the yield information for each reaction.
    Args:
        df (pandas.DataFrame): DataFrame containing reaction strings and RDKit reaction objects.
    The DataFrame is expected to have columns 'mols_pro' (lists of RDKit molecule objects for products) and 'rxn_str' (reaction strings).
    The function identifies biphenol products using a SMARTS pattern and extracts their SMILES representations.
    It also retrieves the yield information from the reaction strings.
    The function adds two new columns to the DataFrame:
        - 'smiles_pro': SMILES strings of the products.
        - 'ylds': Yields of the reactions.
    Returns:
        pandas.DataFrame: DataFrame with additional columns 'smiles_pro' and 'ylds'.
    """
    biphenol_smarts = '[OH]c1c(-c2c([OH])cccc2)cccc1'
    biphenol_mol = Chem.MolFromSmarts(biphenol_smarts)
    smiles_pros = []
    yld = []
    for i in range(len(df)):
        try:
            for j in range(len(df['mols_pro'][i])):
                mol = Chem.MolFromSmiles(Chem.MolToSmiles(df['mols_pro'][i][j]))
                if mol.HasSubstructMatch(biphenol_mol):
                    smiles_pro = Chem.MolToSmiles(mol)
                    smiles_pros.append(smiles_pro)
                    break
        except:
            smiles_pros.append('error')
        
        if i+1 != len(smiles_pros):
            smiles_pros.append('error')
        
        rxn_n = df['rxn_str'][i].split('\n')
        for k in [str(j+1)]:
            try:
                yld.append(int(rxn_n[rxn_n.index('$DTYPE RXN:VAR(1):PRO('+k+'):YIELD')+1].strip('$DATUM ')))
            except:
                yld.append(np.nan)
    
    df['smiles_pro'] = smiles_pros
    df['ylds'] = yld
    return df

def get_reference(df):
    """
    Extracts references title from the reaction strings in the DataFrame.
    This function processes the reaction strings in the DataFrame to extract title of references which reported each reaction.
    Args:
        df (pandas.DataFrame): DataFrame containing reaction strings.
    The DataFrame is expected to have a column 'rxn_str' containing reaction strings in a specific format.
    The function searches for the reference title in the reaction strings and extracts it.
    If the reference title is not found, it appends 'error' to the list.
    Returns:
        pandas.DataFrame: DataFrame with an additional column 'reference' containing the extracted reference titles.
    The function handles cases where the extraction fails by appending 'error' to the list.
    """
    references=[]
    for i in range(len(df)):
        rxn_n=df['rxn_str'][i].split('\n')
        try:
            j=1
            t1=str(rxn_n[rxn_n.index("$DTYPE RXN:VAR(1):REFERENCE(1):TITLE")+j].replace('$DATUM ',''))
            j+=1
            while "$DTYPE" not in str(rxn_n[rxn_n.index("$DTYPE RXN:VAR(1):REFERENCE(1):TITLE")+j]):
                tj=str(rxn_n[rxn_n.index("$DTYPE RXN:VAR(1):REFERENCE(1):TITLE")+j])
                t1=t1+' '+tj
                j+=1
            references.append(t1)
        except:
            references.append('error')
    df['reference']=references
    return df

def main(rdf_folder, result_folder):
    """
    Args:
        rdf_folder (str): Name of the folder containing the RDF file.
        result_folder (str): Name of the folder where the result Excel file will be saved.
    """
    current_file = Path(__file__).resolve()
    parent_dir = current_file.parent.parent
    rdf_path = parent_dir / f'{rdf_folder}'
    result_path = parent_dir / f'{result_folder}'
    rdf_list = glob.glob(f"{rdf_path}/*.rdf")
    df = read_rdf(rdf_list)
    df = rxn_analyze(df)
    df = get_conditions(df)
    df = get_smiles_rea(df)
    df = get_smiles_pro_yld(df)
    df = get_reference(df)
    df = df[df['smiles_pro']!='error']
    df.to_excel(f'{result_path}/after_rdf_analyze.xlsx', index=False)

if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))
    with open("rdf_analyze.yaml", mode="rb") as yml:
        settings = yaml.safe_load(yml)
    
    rdf_folder = settings["rdf_folder"]
    result_folder = settings["result_folder"]

    main(rdf_folder, result_folder)