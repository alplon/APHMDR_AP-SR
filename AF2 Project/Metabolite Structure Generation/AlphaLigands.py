import os
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.rdmolfiles import MolToMolFile
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt

file_in = Path('')#inputpath
dir_out = Path('')#outputpath

df = pd.read_csv(file_in)

os.chdir(dir_out)
print('length is: ', len(df))

with open("test.txt", "a+") as met:
    met.write('YMDB, MOLECULE_NAME, AVERAGE_MASS, SMILES\n')

for counter, molecule in enumerate(df['SMILES']):
    mol = Chem.MolFromSmiles(molecule)
    HMDB = df['MET_ID'][counter]
    Name = df['NAME'][counter]
    Smiles = df['SMILES'][counter]
    try:
        flat = Compute2DCoords(mol)
    except:
        print('this molecule had an error', HMDB,Name)
        continue
    smiles_to_mass = Chem.MolFromSmiles(Smiles)
    molar_mass = Descriptors.MolWt(smiles_to_mass)
        

    if 300 < molar_mass < 800:
        average_mass = str(molar_mass)
        mol_H = Chem.AddHs(mol)
        mol_H_3D = AllChem.EmbedMolecule(mol_H)
        mol_ReH = Chem.RemoveHs(mol_H)
        mol_ReH.SetProp('_Name', Name)
        filename = (HMDB + '.mol')
        MolToMolFile(mol_H,filename)
        with open("testligand.txt", "a+") as met:
            met.write(HMDB + ', ' + Name + ', '+ average_mass + ', ' + Smiles + '\n')
                
            
#this visualization block works in Jupyter, PyCharm not so much
#Draw.MolToMPL(mol2, size=(200, 200))
#plt.show(block=False)
#plt.close()
