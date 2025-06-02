# pulearn_phenol

## Dataset_construction
In order to create a list of **Synthesized Group (SG)** phenols, **rdf_analyze.py** was run on each RDF file in the **RDF_folder**. For each **after_rdf_analyze.xlsx** output by **rdf_analyze.py**, **Preprocessing_SG.py** was executed to output **phenols_SG.xlsx**. The .py file was executed by entering the necessary information in the corresponding .yaml file and then executing the following command.

```bash
python Dataset_construction/rdf_analyze/rdf_analyze.py
python Dataset_construction/Preprocessing_SG/Preprocessing_SG.py
```
python 
## Descriptor_calculation

## Experimental_data
CSV files of the experimental data used in this study were compiled.
## PU_learning
The functions in **Code.py**, **Dimension_Reduction.py** and **result_visualization.py** are used to output the results in **output.ipynb**.
