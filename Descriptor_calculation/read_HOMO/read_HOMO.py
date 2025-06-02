"""Read HOMO orbital energies from the log files of Gaussian single point calculation.

- Convert the unit of HOMO orbital energies from atomic units (a.u.) to eV.
- Read SMILES if that is written in the log file.
- If SMILES could not be read from the log files, read from the gjf files.
- Exclude imaginary frequencies if they are found in the freq calculation.
"""


import glob
import os

from natsort import natsorted
import pandas as pd
import yaml


def is_normal_frequencies(filename):
    """Determine imaginary frequencies.
    
    Args:
        filename (str): Path of the Gaussian log file for freq calculation

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


def main(folder, file1, file2, csv, hartree2ev=27.2114):
    """
    Args:
        folder (str): Path of the folder where the Gaussian log file is saved.
        file1 (str): Name of the Gaussian log file for freq calculation
                     (The numbers were replaced by *).
        file2 (str): Name of the Gaussian log file for single point calculation
                     (The numbers were replaced by *).
        csv (str): Destination path.
        hartree2ev (float, optional): Hartree to eV conversion factor. Defaults to 27.2114.
    """
    
    os.makedirs(os.path.dirname(csv), exist_ok=True)
    freq_files = natsorted(glob.glob(f"{folder}/{file1}"))

    index_smiles_HOMO = []
    unknown_smiles = []
    uncalculated_HOMO = []
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
                        with open(gjf_file, mode="r") as gjf:
                            gjf_lines = gjf.read().splitlines()
                            gjf_smiles_line = [line for line in gjf_lines if "Smiles" in line]
                            smiles = gjf_smiles_line[0].split(":")[1]
                    except:
                        smiles = None
                        unknown_smiles.append(os.path.basename(energy_file))

                if "Normal termination of Gaussian 16" in lines[-1]:
                    HOMO_lines = [line for line in lines if "Alpha  occ. eigenvalues" in line]
                    try:
                        HOMO_energy = float(HOMO_lines[-1].split()[-1]) * hartree2ev
                    except:
                        HOMO_energy = None
                        uncalculated_HOMO.append(os.path.basename(energy_file))
                        print("HOMO was not calculated!")
                else:
                    HOMO_energy = None
                    uncalculated_HOMO.append(os.path.basename(energy_file))
                    print("Error in single point calculation!")

                index_smiles_HOMO.append([index, smiles, HOMO_energy])

        elif is_normal_frequencies(freq_files[i]) == 1:
            imaginary_frequencies.append(os.path.basename(freq_files[i]))
            print("Imaginary frequency!")

        else:
            uncalculated_frequencies.append(os.path.basename(freq_files[i]))
            print("Freq calculation has not been done!")

    df = pd.DataFrame(index_smiles_HOMO, columns=["Index", "SMILES", "HOMO"])
    df.dropna(subset=["HOMO"], inplace=True)
    df.to_csv(csv, index=False)

    txt = csv.replace(".csv", ".txt")
    with open(txt, mode="w") as f:
        f.write(f"Path: {folder}\n\n")
    
        if len(unknown_smiles) != 0:
            f.write("SMILES unknown:\n")
            for name in unknown_smiles:
                f.write(f"{name}\n")
            f.write("\n")
    
        if len(uncalculated_HOMO) != 0:
            f.write("Error:\n")
            for name in uncalculated_HOMO:
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
    with open("read_HOMO_settings.yaml", mode="rb") as yml:
        settings = yaml.safe_load(yml)
    
    folder = settings["folder"]
    file1 = settings["file1"]
    file2 = settings["file2"]
    csv = settings["csv"]

    main(folder, file1, file2, csv)
