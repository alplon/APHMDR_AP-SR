#Author: Alexander P. Plonski 2/15/2022
#This script filters rmsd changes between cross (exp-predicted) and self (exp-exp or predicted-predicted comparisons and writes it into a csv file) based upon rotatable bonds(#line 49). See '/Analysis Output' for example of csv output. 

#Works best in jupyter notebook or a conda rdkit environment
# %%
import os
import glob
import pandas as pd
import csv
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds

# %%
#data directories variables to specify type (useful for when writing into csv and keep track of data type)
root_dir = os.getcwd()
exp_exp = '' #experimental-experimental directory 'Trial Data/Distance Data/experimental-experimental_csv'
af_af = '' #alphafold-alphafold directory 'Trial Data/Distance Data/alphafold-alphafold_csv'
exp_af = '' #experimental-alphafold directory 'Trial Data/Distance Data/experimental-alphafold_Csv'
protein_directory = '' #uses to get protein names, Directory where 'Target - Exp and  Predicted Models.zip' is

t_exp = 'Exp-Exp'
t_af = 'Af-Af'
t_exp_af = 'Exp-Af'

# %%
#grab all distance file for naming purposes later
os.chdir(os.path.join(exp_af, '1'))
file_names = glob.glob('T*')
os.chdir(root_dir)

# %%
#setup csv file which will be used for plotting data
csvcombined = []
csv_file = '/Analysis Output/...'#output csv name

with open(csv_file, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['Protein', 'RMSD [$\AA$]', 'Comparison'])

# %%
def distance_calcandwrite(data1, type):
    for counter, ligand in enumerate(data1['Ligand']):
        with open('/AF2 Project/Metabolite Structure Generation/Selected_ALL_Metabolites.txt', 'r') as file:
            for line in file:
                split_line = line.split()
                met_id = split_line[0]
                met_id = met_id.strip(',')
                smiles = split_line[-1]
                if met_id == ligand:
                    rot_bonds = CalcNumRotatableBonds(Chem.MolFromSmiles(smiles))
                    if rot_bonds == 9:
                        rmsd = data1['RMSD'][counter]
                        print(rmsd)
                        combined_gdt_name = str(gdt_score) + '\n(' + protein_name +')'
                        protein_value = [combined_gdt_name, rmsd, type]
                        csvcombined.append(protein_value)
    with open(csv_file, 'a+') as f:
        writer = csv.writer(f)
        writer.writerows(csvcombined)
    csvcombined.clear()

# %%
#protein normalization based on size, to get rid of normalization, change largest_distance equal to '1'. And run distance_calcandwrite for all trials off al comparisons
#there is incosnistency in naming files for exp-exp and af-af when compared to exp-af, the lines containing exp_name and af_name correct this issue
for fn in file_names:

    protein_name = (fn.split('_'))[0]

    alpha_result_dir = os.path.join (protein_directory, protein_name, 'AlphaFold_Results.txt')#file located in protein_directory for each protein. Contains GDT_TS of model used. 
    with open(alpha_result_dir, 'r') as file:
        for i, line in enumerate(file):
            if i == 1:
                best_model_line = line
                line_split = best_model_line.split()
                gdt_score = line_split[3]


    exp_af_trials = glob.glob(exp_af + '/*')
    for trial in exp_af_trials:
        #print(fn)
        casp_df = pd.read_csv(os.path.join(trial,fn))
        distance_calcandwrite(casp_df, t_exp_af)
        #print(trial)

    trials = ['1-2', '1-3', '2-3']
    for trial in trials:
        exp_name = fn[:-21] + 'Target' + fn[-12:-4] + '_experimental' + trial + fn[-4:]
        #print(exp_name)
        casp_df = pd.read_csv(os.path.join(exp_exp, trial, exp_name))
        distance_calcandwrite(casp_df, t_exp) 

    trials = ['1-2', '1-3', '2-3']
    for trial in trials:
        af_name = fn[:-4] + '_alphafold' + trial + fn[-4:]
        #print(af_name)
        casp_df = pd.read_csv(os.path.join(af_af, trial, af_name))
        distance_calcandwrite(casp_df, t_af)      

# %%



