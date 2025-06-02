"""Read Mulliken charges from the log files of Gaussian single point calculation.

- Read SMILES if that is written in the log file.
- If SMILES could not be read from the log files, read from the gjf files.
- Exclude imaginary frequencies if they are found in the freq calculation.
"""


import glob
import os
from pathlib import Path

from natsort import natsorted
import pandas as pd
from rdkit import Chem
import yaml


def is_normal_frequencies(filename):
    """Determine imaginary frequencies.
    
    Args:
        filename (str): Path of the Gaussian log file for freq calculation.

    Returns:
        int: 0 if there are no imaginary frequencies,
             1 if there are,
             2 if the freq calculation has not been done.
    """
    try:
        with open(filename, mode="r") as f:
            lines = f.readlines()
            if "Normal termination of Gaussian 16" in lines[-1]:
                frequencies_lines = [line for line in lines if "Frequencies" in line]

                if len(frequencies_lines) == 0:
                    return 2

                for l in frequencies_lines:
                    splited = l.split()
                    values = splited[2:]
                    values = [float(v) for v in values]
                    for v in values:
                        if v < 0:
                            return 1
                return 0
            
            else:
                return 2
    
    except:
        return 2


def get_indexes_fromCSV(smiles, df_indexes):
    """Get the index of atoms that compose the phenol skeleton from DataFrame.

    Args:
        smiles (str): SMILES of phenols.
        df_indexes (DataFrame): DataFrame of get_indexes_***.csv by tools.py

    Returns:
        list: [O, next_O, o_1, m_1, p, m_2, o_2]
    """
    smiles_index = df_indexes.index[df_indexes["smiles"] == smiles].tolist()
    smiles_index = smiles_index[0]
    
    return [df_indexes["O"][smiles_index],
            df_indexes["next_O"][smiles_index],
            df_indexes["o_1"][smiles_index],
            df_indexes["m_1"][smiles_index],
            df_indexes["p"][smiles_index],
            df_indexes["m_2"][smiles_index],
            df_indexes["o_2"][smiles_index]]


def main(folder, file1, file2, file3, csv):
    """
    Args:
        folder (str): Path of the folder where the Gaussian log file is saved.
        file1 (str): Name of the Gaussian log file for freq calculation
                     (The numbers were replaced by *).
        file2 (str): Name of the Gaussian log file for single point calculation
                     (The numbers were replaced by *).
        file3 (str): Name of the get_indexes_***.csv by tools.py
        csv (str): Destination path.
    """

    current_file = Path(__file__).resolve()
    parent_dir = current_file.parent.parent
    indexes_CSV_path = parent_dir/"tools"/f"{file3}"
    df_indexes = pd.read_csv(indexes_CSV_path)

    os.makedirs(os.path.dirname(csv), exist_ok=True)
    freq_files = natsorted(glob.glob(f"{folder}/{file1}"))
    
    index_smiles_charges = []
    unknown_smiles = []
    unknown_charges = []
    imaginary_frequencies = []
    uncalculated_frequencies = []

    for i in range(len(freq_files)):
        print(freq_files[i])

        if is_normal_frequencies(freq_files[i]) == 0:
            index = os.path.basename(freq_files[i]).split("_")[0]
            energy_file = f"{folder}/{file2.replace('*', index)}"
            print(energy_file)

            with open(energy_file, mode="r") as f:
                lines = f.read().splitlines()

                smiles_lines_index = []
                smiles_lines = []
                for j, line in enumerate(lines):
                    if "Smiles" in line:
                        smiles_lines_index.append(j)
                try:
                    k = 0
                    while True:
                        if "---" not in lines[smiles_lines_index[0]+k]:
                            smiles_lines.append(lines[smiles_lines_index[0]+k])
                            k += 1
                        else:
                            break
                    smiles = "".join(smiles_lines).replace(" ", "").split(":")[1]
                    if smiles == "None":
                        smiles = None
                        unknown_smiles.append(os.path.basename(energy_file))
                except:
                    try:
                        gjf_file = f"{energy_file.replace('.log', '.gjf')}"
                        with open(gjf_file) as gjf:
                            gjf_lines = gjf.read().splitlines()
                            gjf_smiles_line = [line for line in gjf_lines if "Smiles" in line]
                            smiles = gjf_smiles_line[0].split(":")[1]
                    except: 
                        smiles = None
                        unknown_smiles.append(os.path.basename(energy_file))

                mol = Chem.MolFromSmiles(smiles)
                smiles_checked = Chem.MolToSmiles(mol)

                if "Normal termination of Gaussian 16" in lines[-1]:
                    charges_lines_index = []
                    charges_lines = []
                    for j, line in enumerate(lines):
                        if "Mulliken charges" in line:
                            charges_lines_index.append(j)
                    try:
                        k = 2
                        while True:
                            if "Sum of Mulliken charges" not in lines[charges_lines_index[0]+k]:
                                charges_lines.append(lines[charges_lines_index[0]+k])
                                k += 1
                            else:
                                break

                        C_rdkit_indexes = get_indexes_fromCSV(smiles_checked, df_indexes)
                        C_gaussian_indexes = [x + 1 for x in C_rdkit_indexes]

                        charges = []
                        for l in C_gaussian_indexes:
                            for m in charges_lines:
                                if l == int(m.split()[0]):
                                    charges.append(float(m.split()[-1]))
                                    break
                    except:
                        charges = [None, None, None, None, None, None, None]
                        unknown_charges.append(os.path.basename(energy_file))
                        print("Unable to read charges!")
                else:
                    charges = [None, None, None, None, None, None, None]
                    unknown_charges.append(os.path.basename(energy_file))
                    print("Error in single point calculation!")

                li = sum([[index, smiles], charges], [])
                index_smiles_charges.append(li)

        elif is_normal_frequencies(freq_files[i]) == 1:
            imaginary_frequencies.append(os.path.basename(freq_files[i]))
            print("Imaginary frequency!")

        else:
            uncalculated_frequencies.append(os.path.basename(freq_files[i]))
            print("Freq calculation has not been done!")

    df = pd.DataFrame(index_smiles_charges,
                      columns=["Index", "SMILES",
                               "Charge_O", "Charge_next_O",
                               "Charge_o_1", "Charge_m_1",
                               "Charge_p", "Charge_m_2", "Charge_o_2"])

    df.dropna(subset=["Charge_O", "Charge_next_O",
                      "Charge_o_1", "Charge_m_1",
                      "Charge_p", "Charge_m_2", "Charge_o_2"],
              inplace=True)

    df.to_csv(csv, index=False)

    txt = csv.replace(".csv", ".txt")
    with open(txt, mode="w") as f:
        f.write(f"Path: {folder}\n\n")
    
        if len(unknown_smiles) != 0:
            f.write("SMILES unknown:\n")
            for name in unknown_smiles:
                f.write(f"{name}\n")
            f.write("\n")
    
        if len(unknown_charges) != 0:
            f.write("Error:\n")
            for name in unknown_charges:
                f.write(f"{name}\n")
            f.write("\n")
    
        if len(imaginary_frequencies) != 0:
            f.write("Imaginary frequency:\n")
            for name in imaginary_frequencies:
                f.write(f"{name}\n")
            f.write("\n")
    
        if len(uncalculated_frequencies) != 0:
            f.write("Freq calculation has not been done:\n")
            for name in uncalculated_frequencies:
                f.write(f"{name}\n")

    print(df)

if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))
    with open("read_charges_settings.yaml", mode="rb") as yml:
        settings = yaml.safe_load(yml)
    
    folder = settings["folder"]
    file1 = settings["file1"]
    file2 = settings["file2"]
    file3 = settings["file3"]
    csv = settings["csv"]

    main(folder, file1, file2, file3, csv)
