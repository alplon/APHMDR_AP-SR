#Author: Alexander P. Plonski 2/15/2022
#This script computes distances between cross (exp-predicted) and self (exp-exp or predicted-predicted comparisons and writes it into a csv file.) See '/Analysis Output' for example of csv output. 

# %%
import os
import glob
import pandas as pd
import csv

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
csv_file = '/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/distance_comparisons_all_normalized_trialinfo.csv' #output csv name

with open(csv_file, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['Protein', 'Distance Change [$\AA$]', 'Comparison','Trial'])

# %%
def distance_calcandwrite(data1, type, ctrial):
    delta_distance = abs(data1['ORG Distance'] - data1['MUT Distance'])
    data = delta_distance / largest_distance 

    for value in data:
        combined_gdt_name = str(gdt_score) + '\n(' + protein_name +')'
        protein_value = [combined_gdt_name, value, type, ctrial]
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

    alpha_result_dir = os.path.join (protein_directory, protein_name, 'AlphaFold_Results.txt') #file located in protein_directory for each protein. Contains GDT_TS of model used. 
    with open(alpha_result_dir, 'r') as file:
        for i, line in enumerate(file):
            if i == 1:
                best_model_line = line
                line_split = best_model_line.split()
                gdt_score = line_split[3]


    normalization_file = open("/AF2 Prpject/Analysis Code/Protein Normalization", "r") #Contains greatest distance between two atoms in a protein (essentially 'protein size') for normalization. 
    flag = 0
    index = 0
    
    for line in normalization_file:
        index + 1
        
        if protein_name in line:
            flag = 1
            break
    if flag == 0:
        ignore = 'ignore'
    else:
        largest_distance = 1 #float(line.split(',')[-1]) #*****Uncomment and delete '1' for normalization' Otherwise it will be non-rnormalized distance

    exp_af_trials = glob.glob(exp_af + '/*')
    for trial in exp_af_trials:
        casp_df = pd.read_csv(os.path.join(trial,fn))
        ctrial = t_exp_af + trial[-1]
        distance_calcandwrite(casp_df, t_exp_af,ctrial)
    
    trials = ['1-2', '1-3', '2-3']
    for trial in trials:
        exp_name = fn[:-35] + 'distance' +fn[-16:]
        exp_name = exp_name[:-4] + '_experimental' + trial + exp_name[-4:]
        ctrial = t_exp + trial
        casp_df = pd.read_csv(os.path.join(exp_exp, trial, exp_name))
        distance_calcandwrite(casp_df, t_exp,ctrial)     

    trials = ['1-2', '1-3', '2-3']
    for trial in trials:
        af_name = fn[:-35] + 'distance' +fn[-16:]
        af_name = af_name[:-4] + '_alphafold' + trial + af_name[-4:]
        ctrial = t_af + trial
        casp_df = pd.read_csv(os.path.join(af_af, trial, af_name))
        distance_calcandwrite(casp_df, t_af,ctrial)      



