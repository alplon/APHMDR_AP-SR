# APHMDR_AP-SR
This repository holds data for "Assessing Protein Homology Models with Docking Reproducibility" by Alexander P. Plonski and Scott M. Reed.

Recommended order of viewing project workflow + information about each directory. 

>Metabolite Structure Generation

Includes information about each meatabolite from the Human, Yeast, and E.Coli Metabolomes. Additionally 'AlphaLigands.py' shows how metabolites were selected/generated in 3D structures. Additioanly structure prepation was done with JChem -cxcalc for pH and AutoDockTools for correct file format(.pdbqt). 

>Target - Exp and Predicted Models.zip

Includes information about selected AF2 models from CASP-14 along with Target models (both experimental and predicted)

>Docking Workflow-Output Example

Includes the file organization/workflow of docking with an example contianing a few ligands and one protein and scriptss for file prep, docking, and analysis. 

>Trial Data

Includes all raw distance and rmsd data for cross (experimental-alphafold) and self (experimental-experimental and alphafold-alphafold) comparisons. 

>Analysis Code

Includes distance, lower bulb, and rmsd scripts used to analyze data in '/Trial Data'. Additionally 'data_plot.py' was used to generate figures. 

>Analysis Output

Includes csv output files of processed raw data using scripts in '/Analysis Code'. 
