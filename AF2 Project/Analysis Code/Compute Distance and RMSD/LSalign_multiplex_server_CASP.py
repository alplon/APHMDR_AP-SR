import os
import os.path
import re
import csv
import time

start_time = time.time()

root_dir = os.getcwd()
time_file_path = (root_dir + '/time_file_lsalign.txt')
#print (root_dir)

for fn in os.listdir(root_dir):
    root_file = os.path.join(root_dir, fn)
    #print (root_file)
    if os.path.isdir(root_file):
        os.chdir(root_file)
        native = [fn for fn in os.listdir(".") if re.findall('_Target',fn)]
        mutants = [fn for fn in os.listdir(".") if re.findall('_AlphaFold',fn)]
        for fn in native:
            native_folder = os.path.join(root_file, fn)
            #print (native_folder)
        for fn in os.listdir(native_folder):
            os.chdir(native_folder)
            if os.path.exists("ligandoutput"):
                os.chdir('ligandoutput')
                template_names_mol2 = [fn for fn in os.listdir(".") if fn.endswith('.mol2')]
                os.chdir('..')
        for fn in mutants:
            mutant_folders = os.path.join(root_file, fn)
            mutant_folders_files = os.path.join(mutant_folders, 'LSalignfiles')
            os.chdir(mutant_folders)
            if not os.path.exists('LSalignfiles'):
                os.mkdir('LSalignfiles')
            mutantsx = [fn]
            mutantcsv= fn +"_LSalign.csv"
            if not os.path.isfile(mutantcsv):
                for fn in template_names_mol2:
                    template_name = os.path.splitext(fn)[0]
                    #print(template_names)
                    os.chdir(native_folder)
                    os.chdir('ligandoutput')
                    native_ligandout = os.getcwd()
                    native_mol2 = os.path.join(native_ligandout, fn)
                    print (native_mol2)
                    os.chdir('..')
                    os.chdir('..')
                    os.chdir(mutant_folders)
                    os.chdir('ligandoutput')
                    mutant_ligandout = os.getcwd()
                    mutant_mol2 = os.path.join(mutant_ligandout, fn)
                    print (mutant_mol2)
                    os.chdir('..')
                    os.chdir('LSalignfiles')
                    LScommand = '/opt/lsalign/src/LSalign '+mutant_mol2 +' ' +native_mol2+' >' +' ' +template_name +'.txt'#for Pharmaco
                    #LScommand = '../../../LSalign '+mutant_mol2 +' ' +native_mol2+' >' +' ' +template_name +'.txt' #for mac
                    #print (LScommand)
                    os.system(LScommand)
                    os.chdir('..')
                filenames_nopath = [fn for fn in os.listdir(mutant_folders_files) if fn.endswith(".txt")]
                filenames = []
                for fn in filenames_nopath:
                    filenames.append(os.path.join(mutant_folders_files, fn))
                name_rmsd = []
                def sortligands(filenames):
                    for f in (filenames):
                        outfile = open(f, 'r')
                        data = outfile.readlines()
                        outfile.close()
                        for i, line in enumerate(data):
                            if i == 3:
                                rmsd_line = line
                                numbers = rmsd_line.split()
                                rmsd = str(numbers[5])
                        base = os.path.basename(f)
                        base_name = os.path.splitext(base)[0]
                        #x = base_name.split('_')
                        #z = str(x[1])
                        #if int(x[1]) < 10:
                        #    d = z.rjust(2,'0')
                        #else:
                        #    d = x[1]
                        #ligandname_number = str(x[0]) + '_' + str(d)
                        both = base_name + " " + rmsd
                        name_rmsd.append(both)
                print (mutantsx)
                sortligands(filenames)
                csvcombined = []
                for f in sorted(name_rmsd):
                    sortedRMSD = f
                    ligandcsv, RMSDcsv= sortedRMSD.split()
                    print(ligandcsv + " " + RMSDcsv)
                    combined = [ligandcsv, RMSDcsv]
                    csvcombined.append(combined)
                with open(mutantcsv, 'w') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Ligand', 'RMSD'])
                    writer.writerows(csvcombined)
time_file = open(time_file_path, "w")
time_file.write("--- %s seconds ---" % (time.time() - start_time))
print ('done')
