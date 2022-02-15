#Author: Alexander P. Plonski 2/15/2022
#This script filters distances <=2 Angstroms between cross (exp-predicted) and self (exp-exp or predicted-predicted comparisons and writes it into a csv file.) See '/Analysis Output' for example of csv output. 

#This script was originally run in a jupyter notebok and thus works best. It allows you to run specific blocks. 
# %%
import os
import glob
import pandas as pd
import numpy as np
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
bulb_values = []
def bulb_calcandwrite(data1, type):
    delta_distance = abs(data1['ORG Distance'] - data1['MUT Distance'])

    for value in delta_distance:
        if value <= 2:
            bulb_values.append(value)

    bulb_fraction = len(bulb_values)/1334
    bulb_values.clear()
    return bulb_fraction

# %%
#sort all files by gdt (do this before to make dataframe creation easier for multiple trials)
all_gdt = []
for fn in file_names:
    protein_name = (fn.split('_'))[0]
    #get gdt score
    alpha_result_dir = os.path.join (protein_directory, protein_name, 'AlphaFold_Results.txt')
    with open(alpha_result_dir, 'r') as file:
        for i, line in enumerate(file):
            if i == 1:
                best_model_line = line
                line_split = best_model_line.split()
                gdt_score = line_split[3]
                all_gdt.append(gdt_score)

lists = sorted(zip(*[all_gdt, file_names]))
sorted_gdt, sorted_filenames = list(zip(*lists))

# %%
#Af-Exp
appened_protein_dfs = []
all_bulb_stats = []
for fn in sorted_filenames:
    protein_name = (fn.split('_'))[0]

    appended_dfs = []
    exp_af_trials = glob.glob(exp_af + '/*')
    for trial in exp_af_trials:
        trial_number = trial[-1]
        casp_df = pd.read_csv(os.path.join(trial,fn))
        df1 = pd.DataFrame({'Protein': [protein_name], trial_number:bulb_calcandwrite(casp_df, t_exp_af)},columns=['Protein', '1', '2', '3'])
        appended_dfs.append(df1)
    protein_df = pd.concat(appended_dfs, ignore_index= True)
    protein_df = protein_df.groupby(by=["Protein"], dropna=True).sum()
    appened_protein_dfs.append(protein_df)

    mean = protein_df.mean(axis = 1)
    stdev = protein_df.std(axis = 1)
    bulb_stats = [protein_name, mean.item(), stdev.item()]
    all_bulb_stats.append(bulb_stats)

bulb_stats_df = pd.DataFrame(all_bulb_stats,columns=["Protein","Mean","StDev"]) #output of lower bulb mean and stdev 
bulb_stats_df.to_csv('/path', index=False)

final_df = pd.concat(appened_protein_dfs, ignore_index= False)
final_df.to_csv('/path', index=True) #lowerbulb values for each trial

# %%
#Exp-Exp
appened_protein_dfs = []
all_bulb_stats = []
for fn in sorted_filenames:
    protein_name = (fn.split('_'))[0]

    appended_dfs = []

    trials = ['1-2', '1-3', '2-3']
    for trial in trials:
        exp_name = fn[:-35] + 'distance' +fn[-16:]
        exp_name = exp_name[:-4] + '_experimental' + trial + exp_name[-4:]
        casp_df = pd.read_csv(os.path.join(exp_exp, trial, exp_name))
        df1 = pd.DataFrame({'Protein': [protein_name], trial:bulb_calcandwrite(casp_df, t_exp)},columns=['Protein', '1-2', '1-3', '2-3'])
        appended_dfs.append(df1)
    protein_df = pd.concat(appended_dfs, ignore_index= True)
    protein_df = protein_df.groupby(by=["Protein"], dropna=True).sum()
    appened_protein_dfs.append(protein_df)

    mean = protein_df.mean(axis = 1)
    stdev = protein_df.std(axis = 1)
    bulb_stats = [protein_name, mean.item(), stdev.item()]
    all_bulb_stats.append(bulb_stats)

bulb_stats_df = pd.DataFrame(all_bulb_stats,columns=["Protein","Mean","StDev"])
bulb_stats_df.to_csv('/path', index=False)#output of lower bulb mean and stdev 

final_df = pd.concat(appened_protein_dfs, ignore_index= False)
final_df.to_csv('/path', index=True)#lowerbulb values for each trial

# %%
#Af-Af
appened_protein_dfs = []
all_bulb_stats = []
for fn in sorted_filenames:
    protein_name = (fn.split('_'))[0]

    appended_dfs = []

    trials = ['1-2', '1-3', '2-3']
    for trial in trials:
        af_name = fn[:-35] + 'distance' +fn[-16:]
        af_name = af_name[:-4] + '_alphafold' + trial + af_name[-4:]
        casp_df = pd.read_csv(os.path.join(af_af, trial, af_name))
        df1 = pd.DataFrame({'Protein': [protein_name], trial:bulb_calcandwrite(casp_df,t_af)},columns=['Protein', '1-2', '1-3', '2-3'])
        appended_dfs.append(df1)
    protein_df = pd.concat(appended_dfs, ignore_index= True)
    protein_df = protein_df.groupby(by=["Protein"], dropna=True).sum()
    appened_protein_dfs.append(protein_df)

    mean = protein_df.mean(axis = 1)
    stdev = protein_df.std(axis = 1)
    bulb_stats = [protein_name, mean.item(), stdev.item()]
    all_bulb_stats.append(bulb_stats)

bulb_stats_df = pd.DataFrame(all_bulb_stats,columns=["Protein","Mean","StDev"])
bulb_stats_df.to_csv('/path', index=False)#output of lower bulb mean and stdev 

final_df = pd.concat(appened_protein_dfs, ignore_index= False)
final_df.to_csv('/path', index=True)#lowerbulb values for each trial

# %%
#Af-Exp by GDT
appened_protein_dfs = []
all_bulb_stats = []
for counter, fn in enumerate(sorted_filenames):
    protein_name = (fn.split('_'))[0]
    gdt = sorted_gdt[counter]

    appended_dfs = []
    exp_af_trials = glob.glob(exp_af + '/*')
    for trial in exp_af_trials:
        trial_number = trial[-1]
        casp_df = pd.read_csv(os.path.join(trial,fn))
        df1 = pd.DataFrame({'GDT_TS': [gdt], trial_number:bulb_calcandwrite(casp_df, t_exp_af)},columns=['GDT_TS', '1', '2', '3'])
        appended_dfs.append(df1)
    protein_df = pd.concat(appended_dfs, ignore_index= True)
    protein_df = protein_df.groupby(by=["GDT_TS"], dropna=True).sum()
    appened_protein_dfs.append(protein_df)

    mean = protein_df.mean(axis = 1)
    stdev = protein_df.std(axis = 1)
    bulb_stats = [gdt, mean.item(), stdev.item()]
    all_bulb_stats.append(bulb_stats)

bulb_stats_df = pd.DataFrame(all_bulb_stats,columns=["GDT_TS","Mean","StDev"])
bulb_stats_df.to_csv('/path', index=False)#output of lower bulb mean and stdev


# %%
#Exp-Exp by GDT
appened_protein_dfs = []
all_bulb_stats = []
for counter, fn in enumerate(sorted_filenames):
    protein_name = (fn.split('_'))[0]
    gdt = sorted_gdt[counter]

    appended_dfs = []

    trials = ['1-2', '1-3', '2-3']
    for trial in trials:
        exp_name = fn[:-35] + 'distance' +fn[-16:]
        exp_name = exp_name[:-4] + '_experimental' + trial + exp_name[-4:]
        casp_df = pd.read_csv(os.path.join(exp_exp, trial, exp_name))
        df1 = pd.DataFrame({'GDT_TS': [gdt], trial:bulb_calcandwrite(casp_df, t_exp)},columns=['GDT_TS', '1-2', '1-3', '2-3'])
        appended_dfs.append(df1)
    protein_df = pd.concat(appended_dfs, ignore_index= True)
    protein_df = protein_df.groupby(by=["GDT_TS"], dropna=True).sum()
    appened_protein_dfs.append(gdt)

    mean = protein_df.mean(axis = 1)
    stdev = protein_df.std(axis = 1)
    bulb_stats = [gdt, mean.item(), stdev.item()]
    all_bulb_stats.append(bulb_stats)

bulb_stats_df = pd.DataFrame(all_bulb_stats,columns=["GDT_TS","Mean","StDev"])
bulb_stats_df.to_csv('/path', index=False)#output of lower bulb mean and stdev


# %%
#Af-Af by GDT
appened_protein_dfs = []
all_bulb_stats = []
for counter, fn in enumerate(sorted_filenames):
    protein_name = (fn.split('_'))[0]
    gdt = sorted_gdt[counter]

    appended_dfs = []

    trials = ['1-2', '1-3', '2-3']
    for trial in trials:
        af_name = fn[:-35] + 'distance' +fn[-16:]
        af_name = af_name[:-4] + '_alphafold' + trial + af_name[-4:]
        casp_df = pd.read_csv(os.path.join(af_af, trial, af_name))
        df1 = pd.DataFrame({'GDT_TS': [gdt], trial:bulb_calcandwrite(casp_df,t_af)},columns=['GDT_TS', '1-2', '1-3', '2-3'])
        appended_dfs.append(df1)
    protein_df = pd.concat(appended_dfs, ignore_index= True)
    protein_df = protein_df.groupby(by=["GDT_TS"], dropna=True).sum()
    appened_protein_dfs.append(gdt)

    mean = protein_df.mean(axis = 1)
    stdev = protein_df.std(axis = 1)
    bulb_stats = [gdt, mean.item(), stdev.item()]
    all_bulb_stats.append(bulb_stats)

bulb_stats_df = pd.DataFrame(all_bulb_stats,columns=["GDT_TS","Mean","StDev"])
bulb_stats_df.to_csv('/path', index=False)#output of lower bulb mean and stdev
