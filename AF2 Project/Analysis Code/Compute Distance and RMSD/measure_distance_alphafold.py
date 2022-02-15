import os
import os.path
import re
import csv
import glob
from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import pandas as pd
from numpy import mean

root_dir = os.getcwd()
pmol = PandasMol2()
ppdb = PandasPdb()

trial_1 = 'casp2'
trial_2 = 'casp3'
trial_3 = 'casp4'

csv_location = '/home/alexanderplonski/alphafold-alphafold_csv'
csv1_2 = '1-2'
csv1_3 = '1-3'
csv2_3 = '2-3'


def measure_distance(data1, data2, csvtype):
    casp_trials = [data1, data2]
    data1_directory = os.path.join(root_dir, data1)
    data2_directory = os.path.join(root_dir, data2)
    
    for dir in os.listdir(data1_directory):
        target_directory = os.path.join(data1_directory, dir)
        
        if os.path.isdir(target_directory):
            target_name = dir
            os.chdir(target_directory)
            experimental1_directory = os.path.join(target_directory, dir + '_AlphaFold')
            experimental2_directory = os.path.join(data2_directory, dir, dir + '_AlphaFold')

            os.chdir(experimental1_directory)
            
            if os.path.exists("ligandoutput"):
                os.chdir('ligandoutput')
                template_names_mol2 = glob.glob('*.mol2')
                os.chdir('..')
            csv = csv_location + '/' + csvtype + '/' + target_name + '_distance_alphacarbon_alphafold' + csvtype + '.csv'
            experimental_pdb = (dir + '_AlphaFold.pdb')
            experimental1_pdb_path = os.path.join(experimental1_directory, experimental_pdb)
            experimental2_pdb_path = os.path.join(experimental2_directory, experimental_pdb)
            
            #print(experimental1_pdb_path)

            target_appended_smalldistance= []
            appended_CAdistance = []
            if not os.path.isfile(csv):
                for fn in template_names_mol2:
                    template_name = os.path.splitext(fn)[0]
                    experimental_mol2path = os.path.join(experimental1_directory, 'ligandoutput', fn)
                    os.chdir(experimental1_directory)
                    os.chdir('ligandoutput')
                    with open(template_name+'_out.pdbqt',"r") as outfile:
                        data = outfile.readlines()
                    for line in data:
                        if 'REMARK VINA RESULT' in line:
                            score_line = line
                            best_score = score_line.split()[3]
                            #print(best_score)
                            break
                    
                    
                    ligandmol2 = pmol.read_mol2(experimental_mol2path)
                    ATOM = pmol.df

                    x_coords = ATOM['x']
                    y_coords = ATOM['y']
                    z_coords = ATOM['z']

                    x_center = mean(x_coords)
                    y_center = mean(y_coords)
                    z_center = mean(z_coords)
                    
                    native_ligand_center = x_center,y_center,z_center
                
                
                    proteinpdb = ppdb.read_pdb(experimental1_pdb_path)
                    ATOM = ppdb.df['ATOM']=ppdb.df['ATOM']

                    
                    
                    target_all_distances = []
                    target_residue_names = []
                    target_residue_numbers = []
                    for counter, atom in enumerate(ATOM['atom_name']):
                        if atom == 'CA':
                        
                            target_residue_names.append(ATOM['residue_name'][counter])
                            target_residue_numbers.append(ATOM['residue_number'][counter])
                        
                            CA_x = ATOM['x_coord'][counter]
                            CA_y = ATOM['y_coord'][counter]
                            CA_z = ATOM['z_coord'][counter]
                            
                            compute_distance = ((CA_x - x_center)**2 + (CA_y - y_center)**2 + (CA_z - z_center)**2)**0.5
                            target_all_distances.append(compute_distance)
                            
                    df1 = pd.DataFrame({'ORG Residue': target_residue_names, 'ORG Residue Number': target_residue_numbers, 'ORG Residue Atom': 'CA', 'ORG Distance': target_all_distances, 'ORG Ligand': template_name, 'ORG Protein': target_name+"_AlphaFold_" + data1, 'ORG Score': best_score, 'ORG Drug X': x_center, 'ORG Drug Y': y_center, 'ORG Drug Z': z_center}, columns=['ORG Protein','ORG Ligand','ORG Residue','ORG Residue Number', 'ORG Residue Atom', 'ORG Distance', 'ORG Score', 'ORG Drug X', 'ORG Drug Y', 'ORG Drug Z'])
                    df1_smallest_distance = df1.nsmallest(1, 'ORG Distance')
                    reference_residue_name = df1_smallest_distance['ORG Residue'].to_numpy()
                    reference_residue_number = df1_smallest_distance['ORG Residue Number'].to_numpy()
                    target_appended_smalldistance.append(df1_smallest_distance)
             
             
                    #proteinpdb = ppdb.read_pdb(experimental1_pdb_path)

                    template_name = os.path.splitext(fn)[0]
                    experimental_mol2path = os.path.join(experimental2_directory, 'ligandoutput', fn)
                    os.chdir(experimental2_directory)
                    os.chdir('ligandoutput')
                    with open(template_name+'_out.pdbqt',"r") as outfile:
                        data = outfile.readlines()
                    for line in data:
                        if 'REMARK VINA RESULT' in line:
                            score_line = line
                            best_score = score_line.split()[3]
                            #print(best_score)
                            break
                            
                    print(experimental_mol2path)
                    ligandmol2 = pmol.read_mol2(experimental_mol2path)
                    ATOM = pmol.df

                    x_coords = ATOM['x']
                    y_coords = ATOM['y']
                    z_coords = ATOM['z']

                    x_center = mean(x_coords)
                    y_center = mean(y_coords)
                    z_center = mean(z_coords)

                    mutant_ligand_center = x_center,y_center,z_center
                    
                    proteinpdb = ppdb.read_pdb(experimental2_pdb_path)
                    ATOM = ppdb.df['ATOM']=ppdb.df['ATOM']
 
                    
                    alpha_all_distances = []
                    alpha_atom_names = []
                    alpha_residue_names = []
                    alpha_residue_numbers = []
                    for counter, number in enumerate(ATOM['residue_number']):
                        if number == reference_residue_number:
                            
                            alpha_atom_names.append(ATOM['atom_name'][counter])
                            alpha_residue_names.append(ATOM['residue_name'][counter])
                            alpha_residue_numbers.append(ATOM['residue_number'][counter])
                        
                            CA_x = ATOM['x_coord'][counter]
                            CA_y = ATOM['y_coord'][counter]
                            CA_z = ATOM['z_coord'][counter]
                            
                            compute_distance = ((CA_x - x_center)**2 + (CA_y - y_center)**2 + (CA_z - z_center)**2)**0.5
                            alpha_all_distances.append(compute_distance)

                    df1 = pd.DataFrame({'MUT Residue': alpha_residue_names, 'MUT Residue Number': alpha_residue_numbers, 'MUT Residue Atom': alpha_atom_names, 'MUT Distance': alpha_all_distances, 'MUT Ligand': template_name, 'MUT Protein': target_name+"_AlphaFold_" + data2, 'MUT Score': best_score, 'MUT Drug X': x_center, 'MUT Drug Y': y_center, 'MUT Drug Z': z_center}, columns=['MUT Protein','MUT Ligand','MUT Residue','MUT Residue Number', 'MUT Residue Atom', 'MUT Distance', 'MUT Score', 'MUT Drug X', 'MUT Drug Y', 'MUT Drug Z'])
                    df1_CA_atom = df1.loc[df1['MUT Residue Atom'] == 'CA']
                    appended_CAdistance.append(df1_CA_atom)
                    
                target_smalldistance_df = pd.concat(target_appended_smalldistance,ignore_index=True)
                CAdistance_df = pd.concat(appended_CAdistance,ignore_index=True)
                
                target_alpha_df= pd.concat([target_smalldistance_df, CAdistance_df], sort=False, axis=1,)
                target_alpha_sorted_df = target_alpha_df.sort_values(by=['ORG Ligand'], ascending=True)

                os.chdir(root_dir)
                target_alpha_dftocsv = target_alpha_sorted_df.to_csv(csv, index=False)


measure_distance(trial_1, trial_2, csv1_2)
measure_distance(trial_1, trial_3, csv1_3)
measure_distance(trial_2, trial_3, csv2_3)


