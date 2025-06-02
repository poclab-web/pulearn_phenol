# pulearn_phenol

## Dataset_construction
In order to create a list of **Synthesized Group (SG)** phenols, **rdf_analyze.py** was run on each RDF file in the **RDF_folder**. For **after_rdf_analyze.xlsx** output by **rdf_analyze.py**, **Preprocessing_SG.py** was executed to output **phenols_SG.xlsx**. The .py file was executed by entering the necessary information in the corresponding .yaml file and then executing the following command.

```bash
python Dataset_construction/rdf_analyze/rdf_analyze.py
python Dataset_construction/Preprocessing_SG/Preprocessing_SG.py
```

In order to create a list of **Biological Group (BG)** phenols, **Preprocessing_BG.py** was run on each SDF file in the **SDF_BG**.  Results were summarized in **phenol_BG.xlsx** . The **Preprocessing_BG.py** file was executed by entering the necessary information in **Preprocessing_BG.yaml** and then executing the following command.

```bash
python Dataset_construction/Preprocessing_BG/Preprocessing_BG.py
```

In order to create a list of **Reagents Group (RG)** phenols, **Preprocessing_RG.py** was run on each SDF file in the **SDF_RG**.  Results were summarized in **phenol_RG.csv** . The **Preprocessing_RG.py** file was executed by entering the necessary information in **Preprocessing_RG.yaml** and then executing the following command. **SDF_RG** could not be uploaded to this repository because the folder was too large.
```bash
python Dataset_construction/Preprocessing_RG/Preprocessing_RG.py
```
## Descriptor_calculation
MolLogP, MolWt, and Number of heteroatoms were calculated by executing the following commands.

```bash
python Descriptor_calculation/rdkit_descriptors/calc_rdkit_descriptors.py
```

The HOMO and Mulliken charge values were extracted from the Gaussian.logfile by entering the necessary information in the corresponding .yaml files and executing the following commands, respectively.

```bash
python Descriptor_calculation/read_HOMO/read_HOMO.py
python Descriptor_calculation/read_charges/read_charges.py
```

The Sterimol values (L, B1, B5) for the substituents were calculated by executing the following commands
```bash
python Descriptor_calculation/sterimol/calculate_sterimol.py {Path of the folder where the result files (Gaussian.logfile) of structural optimization of substituents are stored.} {Path to CSV file of calculation results}
```
Mulliken chargeを

## Experimental_data
CSV files of the experimental data used in this study were compiled.
## PU_learning
The functions in **Code.py**, **Dimension_Reduction.py** and **result_visualization.py** are used to output the results in **output.ipynb**.
