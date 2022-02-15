import os
import os.path
import re
import csv
from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import pandas as pd
from numpy import mean

root_dir = os.getcwd()
pmol = PandasMol2()
ppdb = PandasPdb()


for fn in os.listdir(root_dir):
    root_file = os.path.join(root_dir, fn)
    if os.path.isdir(root_file):
        os.chdir(root_file)
        native = [fn for fn in os.listdir(".") if re.findall('_Target',fn)]
        mutants = [fn for fn in os.listdir(".") if re.findall('_AlphaFold',fn)]
        for fn in native:
            native_folder = os.path.join(root_file, fn)
        for fn in os.listdir(native_folder):
            os.chdir(native_folder)
            if os.path.exists("ligandoutput"):
                os.chdir('ligandoutput')
                template_names_mol2 = [fn for fn in os.listdir(".") if fn.endswith('.mol2')]
                os.chdir('..')
        for fn in mutants:
            mutant_folders = os.path.join(root_file, fn)
            mutantsx = [fn]
            mutantcsv = fn+"_distscore_alphacarbon.csv"
            proteinandresidue = fn.split('_')
            residue_sitesplit= str(proteinandresidue[1])
            protein_name= str(proteinandresidue[0])
            #print(protein_name)
            mutantpdb = (fn+".pdb")
            nativepdb = (protein_name+"_Target.pdb")
            native_pdbpath = os.path.join(native_folder, nativepdb)
            print(native_pdbpath)
            mutant_pdbpath = os.path.join(mutant_folders, mutantpdb)
            print(mutant_pdbpath)
            os.chdir(mutant_folders)
  
  
            target_appended_smalldistance= []
            appended_CAdistance = []
            if not os.path.isfile(mutantcsv):
                for fn in template_names_mol2:
                    template_name = os.path.splitext(fn)[0]
                    native_mol2path = os.path.join(native_folder, 'ligandoutput', fn)
                    os.chdir(native_folder)
                    os.chdir('ligandoutput')
                    with open(template_name+'_out.pdbqt',"r") as outfile:
                        data = outfile.readlines()
                    for line in data:
                        if 'REMARK VINA RESULT' in line:
                            score_line = line
                            best_score = score_line.split()[3]
                            #print(best_score)
                            break
                    
                    
                    ligandmol2 = pmol.read_mol2(native_mol2path)
                    ATOM = pmol.df

                    x_coords = ATOM['x']
                    y_coords = ATOM['y']
                    z_coords = ATOM['z']

                    x_center = mean(x_coords)
                    y_center = mean(y_coords)
                    z_center = mean(z_coords)
                    
                    native_ligand_center = x_center,y_center,z_center
                
                
                    proteinpdb = ppdb.read_pdb(native_pdbpath)
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
                            
                    df1 = pd.DataFrame({'ORG Residue': target_residue_names, 'ORG Residue Number': target_residue_numbers, 'ORG Residue Atom': 'CA', 'ORG Distance': target_all_distances, 'ORG Ligand': template_name, 'ORG Protein': protein_name+"_Target", 'ORG Score': best_score, 'ORG Drug X': x_center, 'ORG Drug Y': y_center, 'ORG Drug Z': z_center}, columns=['ORG Protein','ORG Ligand','ORG Residue','ORG Residue Number', 'ORG Residue Atom', 'ORG Distance', 'ORG Score', 'ORG Drug X', 'ORG Drug Y', 'ORG Drug Z'])
                    df1_smallest_distance = df1.nsmallest(1, 'ORG Distance')
                    reference_residue_name = df1_smallest_distance['ORG Residue'].to_numpy()
                    reference_residue_number = df1_smallest_distance['ORG Residue Number'].to_numpy()
                    target_appended_smalldistance.append(df1_smallest_distance)


                    proteinpdb = ppdb.read_pdb(mutant_pdbpath)

                    template_name = os.path.splitext(fn)[0]
                    mutant_mol2path = os.path.join(mutant_folders, 'ligandoutput', fn)
                    os.chdir(mutant_folders)
                    os.chdir('ligandoutput')
                    with open(template_name+'_out.pdbqt',"r") as outfile:
                        data = outfile.readlines()
                    for line in data:
                        if 'REMARK VINA RESULT' in line:
                            score_line = line
                            best_score = score_line.split()[3]
                            #print(best_score)
                            break
                            
                    print(mutant_mol2path)        
                    ligandmol2 = pmol.read_mol2(mutant_mol2path)
                    ATOM = pmol.df

                    x_coords = ATOM['x']
                    y_coords = ATOM['y']
                    z_coords = ATOM['z']

                    x_center = mean(x_coords)
                    y_center = mean(y_coords)
                    z_center = mean(z_coords)

                    mutant_ligand_center = x_center,y_center,z_center
                    
                    proteinpdb = ppdb.read_pdb(mutant_pdbpath)
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

                    df1 = pd.DataFrame({'MUT Residue': alpha_residue_names, 'MUT Residue Number': alpha_residue_numbers, 'MUT Residue Atom': alpha_atom_names, 'MUT Distance': alpha_all_distances, 'MUT Ligand': template_name, 'MUT Protein': protein_name+"_AlphaFold", 'MUT Score': best_score, 'MUT Drug X': x_center, 'MUT Drug Y': y_center, 'MUT Drug Z': z_center}, columns=['MUT Protein','MUT Ligand','MUT Residue','MUT Residue Number', 'MUT Residue Atom', 'MUT Distance', 'MUT Score', 'MUT Drug X', 'MUT Drug Y', 'MUT Drug Z'])
                    df1_CA_atom = df1.loc[df1['MUT Residue Atom'] == 'CA']
                    appended_CAdistance.append(df1_CA_atom)
                    
                target_smalldistance_df = pd.concat(target_appended_smalldistance,ignore_index=True)
                CAdistance_df = pd.concat(appended_CAdistance,ignore_index=True)
                
                target_alpha_df= pd.concat([target_smalldistance_df, CAdistance_df], sort=False, axis=1,)
                target_alpha_sorted_df = target_alpha_df.sort_values(by=['ORG Ligand'], ascending=True)
                
                os.chdir(mutant_folders)
                target_alpha_dftocsv = target_alpha_sorted_df.to_csv(mutantcsv, index=False)
                os.chdir(root_dir)
            
print("done")
