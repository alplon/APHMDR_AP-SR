#This code will multiplex through all proteins and create a pdbqt file for them if necessary and run the qvina docking process if necessary. Run in root directory. Last Updated: 01-05-2021 by AP
import os
import os.path
import re
import time

start_time = time.time()

os.system('export PATH=/home/boss/qvina:$PATH')
os.system('export LD_LIBRARY_PATH=/home/boss/qvina/lib')

#9-17 appends all protein folders (mutant and natives) into the variable 'protein'
root_dir = os.getcwd()
time_file_path = (root_dir + '/time_file.txt')
protein = []
for fn in os.listdir(root_dir):
    root_file = os.path.join(root_dir, fn)
    if os.path.isdir(root_file):
        os.chdir(root_file)
        sub_dir = [fn for fn in os.listdir(".") if os.path.isdir(fn)]
        for fn in sub_dir:
            protein.append(fn)

#2-21 gets all mol2 files under the variable 'ligand_names'
os.chdir('..')
ligand_names = [fn for fn in os.listdir(".") if fn.endswith('.mol2')]

#lines 25-30 setup variables for the future lines. lines 31-35 will search for pdbqt files, if hey don't exit, then the ADT command (lines 35,36) will generate a pdbqt from an existing pdb. lines 37-39 check for a ligandoutput folder, if it doesn't exist it will generate one for the protein. lines 41-55 will run the qvina docking command (lines 54,55) against all ligands that are apart of the variable 'ligand_names' (lines 21).
for fp in protein:
    protein_str = str(protein)
    protein_split = fp.split('_')
    protein_base_name = str(protein_split[0])
    protein_path  = os.path.join(root_dir, protein_base_name, fp)
    protein_pdb = (fp + '.pdb')
    os.chdir(protein_path)
    if os.path.isfile(fp + '.pdbqt'):
        print('pdbqt file for ' + fp + 'already exists')
    else:
        print('Converting ' + protein_pdb + ' to pdbqt format')
        pdbqt_command = ('/opt/mgltools/mgltools_x86_64Linux2_1.5.7/bin/pythonsh /opt/mgltools/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r' + protein_pdb)
        os.system(pdbqt_command)
    if not os.path.isdir('ligandoutput'):
        os.mkdir('ligandoutput')
        
    for fn in ligand_names:
        ligand_path = os.path.join(root_dir, fn)
        ligand_splitmol2 = os.path.splitext(fn)[0]
        ligand_pdbqt = (ligand_splitmol2 + '.pdbqt')
        ligand_pdbqt_path = os.path.join(root_dir, ligand_pdbqt)
        ligand_out_pdbqt =(ligand_splitmol2 + '_out.pdbqt')
        ligandoutput_path = os.path.join(protein_path, 'ligandoutput')
        os.chdir(ligandoutput_path)
        #os.system('pwd')
        if not os.path.isfile(ligand_out_pdbqt):
            os.chdir('..')
           # print('docking ' + ligand_splitmol2 + ' with ' + fp)
            os.system('pwd')
            qvina_command = ('qvina-w --config conf.txt --ligand '+ ligand_pdbqt_path + ' --out ligandoutput/' + ligand_out_pdbqt)
            os.system(qvina_command)
        

os.chdir(root_dir)
time_file = open(time_file_path, "w")
time_file.write("--- %s seconds ---" % (time.time() - start_time))
print('done')
