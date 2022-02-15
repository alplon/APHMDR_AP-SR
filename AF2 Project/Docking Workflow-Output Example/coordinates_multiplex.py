import os
import os.path
import re
from biopandas.pdb import PandasPdb
import pandas as pd
import numpy as np
import time

start_time = time.time()

ppdb = PandasPdb()

root_dir = os.getcwd()
protein = []
for fn in os.listdir(root_dir):
    root_file = os.path.join(root_dir, fn)
    if os.path.isdir(root_file):
        os.chdir(root_file)
        sub_dir = [fn for fn in os.listdir(".") if os.path.isdir(fn)]
        for fn in sub_dir:
            protein.append(fn)
            
for fp in protein:
    protein_str = str(protein)
    protein_split = fp.split('_')
    protein_base_name = str(protein_split[0])
    protein_path  = os.path.join(root_dir, protein_base_name, fp)
    protein_pdb = (fp + '.pdb')
    protein_pdbqt = (fp + '.pdbqt')
    protein_path_pdb = os.path.join(protein_path, protein_pdb)
    protein_conf_path = (protein_path + '/conf.txt')
    
    proteinpdb = ppdb.read_pdb(protein_path_pdb)
    ATOM = ppdb.df['ATOM']=ppdb.df['ATOM']


    ref = (0, 0, 0)
    distances = proteinpdb.distance(xyz=ref, records=('ATOM',))
    distances_df = pd.DataFrame({'Distance':distances, 'X':ATOM['x_coord'],'Y':ATOM['y_coord'],'Z':ATOM['z_coord']})


    largestdistance_df = distances_df.nlargest(1, 'Distance')
    x_large_filtered_df = pd.DataFrame(largestdistance_df, columns =['X'])
    y_large_filtered_df = pd.DataFrame(largestdistance_df, columns =['Y'])
    z_large_filtered_df = pd.DataFrame(largestdistance_df, columns =['Z'])
    x_large_array = x_large_filtered_df.to_numpy()
    y_large_array = y_large_filtered_df.to_numpy()
    z_large_array = z_large_filtered_df.to_numpy()


    smallestdistance_df = distances_df.nsmallest(1, 'Distance')
    x_small_filtered_df = pd.DataFrame(smallestdistance_df, columns =['X'])
    y_small_filtered_df = pd.DataFrame(smallestdistance_df, columns =['Y'])
    z_small_filtered_df = pd.DataFrame(smallestdistance_df, columns =['Z'])
    x_small_array = x_small_filtered_df.to_numpy()
    y_small_array = y_small_filtered_df.to_numpy()
    z_small_array = z_small_filtered_df.to_numpy()
    

    x_midpoint = x_small_array + ((x_large_array - x_small_array)/2)
    #print(x_midpoint)
    y_midpoint = y_small_array + ((y_large_array - y_small_array)/2)
    #print(y_midpoint)
    z_midpoint = z_small_array + ((z_large_array - z_small_array)/2)
    #print(z_midpoint)

    x_coords = ATOM['x_coord'].to_numpy()
    y_coords = ATOM['y_coord'].to_numpy()
    z_coords = ATOM['z_coord'].to_numpy()


    x_coord_differences = []
    for item in x_coords:
        x_coord_differences.append(item - x_midpoint)
    x_abs_max = max(x_coord_differences, key = abs)
    x_abs_max_coord = x_midpoint + x_abs_max

    y_coord_differences = []
    for item in y_coords:
        y_coord_differences.append(item - y_midpoint)
    y_abs_max = max(y_coord_differences, key = abs)
    y_abs_max_coord = y_midpoint + y_abs_max

    z_coord_differences = []
    for item in z_coords:
        z_coord_differences.append(item - z_midpoint)
    z_abs_max = max(z_coord_differences, key = abs)
    z_abs_max_coord = z_midpoint + z_abs_max


    x_coord_differences_2 = []
    for item in x_coords:
        x_coord_differences_2.append(item - x_abs_max_coord)
    x_abs_max_2 = max(x_coord_differences_2, key = abs)
    x_abs_max_2_coord = (x_abs_max_coord + x_abs_max_2)
    x_abs_max_list = (x_abs_max_coord, x_abs_max_2_coord)
    x_abs_max_list_sorted = sorted(x_abs_max_list)
    x_upper = x_abs_max_list_sorted[1]
    x_lower = x_abs_max_list_sorted[0]

    y_coord_differences_2 = []
    for item in y_coords:
        y_coord_differences_2.append(item - y_abs_max_coord)
    y_abs_max_2 = max(y_coord_differences_2, key = abs)
    y_abs_max_2_coord = (y_abs_max_coord + y_abs_max_2)
    y_abs_max_list = (y_abs_max_coord, y_abs_max_2_coord)
    y_abs_max_list_sorted = sorted(y_abs_max_list)
    y_upper = y_abs_max_list_sorted[1]
    y_lower = y_abs_max_list_sorted[0]

    z_coord_differences_2 = []
    for item in z_coords:
        z_coord_differences_2.append(item - z_abs_max_coord)
    z_abs_max_2 = max(z_coord_differences_2, key = abs)
    z_abs_max_2_coord = (z_abs_max_coord + z_abs_max_2)
    z_abs_max_list = (z_abs_max_coord, z_abs_max_2_coord)
    z_abs_max_list_sorted = sorted(z_abs_max_list)
    z_upper = z_abs_max_list_sorted[1]
    z_lower = z_abs_max_list_sorted[0]

    x_center_array = x_lower + ((x_upper - x_lower)/2)
    y_center_array = y_lower + ((y_upper - y_lower)/2)
    z_center_array = z_lower + ((z_upper - z_lower)/2)
    
    x_center = x_center_array.item(0)
    y_center = y_center_array.item(0)
    z_center = z_center_array.item(0)

    x_grid_array = x_upper - x_lower + 5
    y_grid_array = y_upper - y_lower + 5
    z_grid_array = z_upper - z_lower + 5
    
    x_grid = x_grid_array.item(0)
    y_grid = y_grid_array.item(0)
    z_grid = z_grid_array.item(0)
    
    
    conf_file = open(protein_conf_path,"w")
    conf_file.write('receptor = ' + protein_pdbqt + "\n")
    conf_file.write("\n")
    conf_file.write('center_x = ' + str(x_center) + "\n")
    conf_file.write('center_y = ' + str(y_center) + "\n")
    conf_file.write('center_z = ' + str(z_center) + "\n")
    conf_file.write('size_x = ' + str(x_grid) + "\n")
    conf_file.write('size_y = ' + str(y_grid) + "\n")
    conf_file.write('size_z = ' + str(z_grid) + "\n")
    conf_file.write("\n")
    conf_file.write('exhaustiveness = 10')

print('done')
print("--- %s seconds ---" % (time.time() - start_time))
