#!/bin/bash

for f in *.mol2; do
b=`basename $f .mol2`
echo Converting mol2 ligands
/opt/mgltools/mgltools_x86_64Linux2_1.5.7/bin/pythonsh /opt/mgltools/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l $f

done
