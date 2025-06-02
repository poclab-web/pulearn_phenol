"""Calculate Sterimol parameters from Gaussian log files.

2H (deuterium) at the end of the substituents.

Outputs the following files
- Calculation results in CSV format (xxxxx.csv)
- Calculation results in text format (xxxxx_log.txt)
- Summary of execution results (xxxxx.txt)
- Input file (xxxxx_input.csv)
"""

import os
import subprocess
import sys

import pandas as pd
from rdkit import Chem


def make_input(txt, csv):
    """Specify the orientation of substituents.

    Args:
        txt (str): Path of the text file in which the SMILES was written (.txt)
        csv (str): Path to save the input file (.csv)
    """
    df = pd.read_table(txt, names=["SMILES"])

    df["ROMol"] = df["SMILES"].map(Chem.MolFromSmiles)
    df["SMILES"] = df["ROMol"].map(Chem.MolToSmiles)
    df["ROMol"] = df["SMILES"].map(Chem.MolFromSmiles)

    li1 = []
    for mol in df["ROMol"]:
        index = mol.GetSubstructMatch(Chem.MolFromSmarts("[2H]*"))
        li1.append(index[0] + 1)
    df.insert(0, "A", li1)

    li2 = []
    for i in range(len(df)):
        a = df.at[i, "ROMol"].GetAtomWithIdx(li1[i] - 1)
        b = a.GetNeighbors()[0].GetIdx()
        li2.append(b + 1)
    df.insert(1, "B", li2)

    df.insert(0, "DFTindex", [i + 1 for i in range(len(df))])
    df.drop(columns=["ROMol"], inplace=True)
    
    print(df)
    df.to_csv(csv, index=False, mode="x")


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


def calculate(folder, csv, txt):
    """Calculate the Sterimol parameters.

    Args:
        folder (str): Path of the folder where the Gaussian log file is saved.
        csv (str): Input file path (.csv)
        txt (str): Path to save the calculation results (.txt)

    Returns:
        tuple: (list, list)
               Names of the file with imaginary frequency,
               Names of the file that has not been freq calculation
    """
    gaussian_log = "*_#_opt_freq_b3lyp_6-31G_d_geom_allcheck_guess_read.log"

    df = pd.read_csv(csv)

    with open(txt, mode="x") as f:
        f.write(f"{txt}\n\n---\n")

    imaginary_frequencies = []
    uncalculated_frequencies = []

    for i in range(len(df)):
        index = df.at[i, "DFTindex"]
        a = df.at[i, "A"]
        b = df.at[i, "B"]
        smiles = df.at[i, "SMILES"]
        filename = gaussian_log.replace("*", str(index))
        file = f"{folder}/{filename}"

        print(index, a, b, filename)
        
        if is_normal_frequencies(file) == 0:
            output = subprocess.run(f"python -m sterimol -a1 {a} -a2 {b} {file}",
                                    shell=True, capture_output=True, text=True)
            with open(txt, mode="a") as f:
                f.write(f"{index}. {smiles}\n")
                f.write(output.stdout)
                f.write("---\n")

        elif is_normal_frequencies(file) == 1:
            imaginary_frequencies.append(filename)
            print("Imaginary frequency!")

        else:
            uncalculated_frequencies.append(filename)
            print("Freq calculation has not been done!")

    return imaginary_frequencies, uncalculated_frequencies


def read(txt, csv_input, csv_output):
    """Read the calculation results.

    Args:
        txt (str): Path of the calculation results in text format (.txt)
        csv_input (str): Input file path (.csv)
        csv_output (str): Destination path (.csv)

    Returns:
        int: Number of substituents for which Sterimol parameters were calculated.
    """
    df_input = pd.read_csv(csv_input)

    index_sterimol = []
    with open(txt, mode="r") as f:
        lines = f.read().splitlines()

        for i in range(len(df_input)):
            index = df_input.at[i, "DFTindex"]
            smiles = df_input.at[i, "SMILES"]

            title_line = None
            sterimol = None
            
            for j, line in enumerate(lines):
                if f"{index}. {smiles}" in line:
                    title_line = j
                    break

            try:
                k = 0
                while True:
                    if "---" not in lines[title_line + k]:
                        k += 1
                        if "Structure" in lines[title_line + k]:
                            sterimol_line = lines[title_line + k + 1]
                            break
                    else:
                        break
                sterimol = sterimol_line.split()
                index_sterimol.append([index, sterimol[1], sterimol[2], sterimol[3]])

            except:
                pass

    df_sterimol = pd.DataFrame(index_sterimol, columns=["DFTindex", "L", "B1", "B5"])
    df_output = pd.merge(df_input, df_sterimol, on="DFTindex", how="inner")

    print(df_output)
    df_output.to_csv(csv_output, index=False, mode="x")

    return len(df_output)


def main(folder, csv):
    """
    Args:
        folder (str): Path of the folder where the Gaussian log file is saved.
        csv (str): Destination path (.csv)
    """
    directory = os.path.dirname(csv)
    os.makedirs(directory, exist_ok=True)

    print("Specifying the orientation of substituents...")
    txt_smiles = f"{folder}/Smiles.txt"
    csv_input = csv.replace(".csv", "_input.csv")
    make_input(txt_smiles, csv_input)

    print("Calculating the Sterimol parameters...")
    txt_log = csv.replace(".csv", "_log.txt")
    imaginary_uncalculated = calculate(folder, csv_input, txt_log)

    print("Writing the calculation results...")
    n = read(txt_log, csv_input, csv)

    txt = csv.replace(".csv", ".txt")
    with open(txt, mode="x") as f:
        f.write(f"Referenced Gaussian output file: {folder}\n")
        f.write(f"Calculation results in CSV format: {csv}\n")
        f.write(f"Number of substituents: {n}\n\n")

        if len(imaginary_uncalculated[0]) != 0:
            f.write("Imaginary frequency:\n")
            for name in imaginary_uncalculated[0]:
                f.write(f"{name}\n")
            f.write("\n")

        if len(imaginary_uncalculated[1]) != 0:
            f.write("Freq calculation has not been done:\n")
            for name in imaginary_uncalculated[1]:
                f.write(f"{name}\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
