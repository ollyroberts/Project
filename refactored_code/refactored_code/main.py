#!/usr/bin/python3
import subprocess
import argparse
import core_functions as core

"""
for one_helix shell 
#!/bin/sh
for file in *.pdb
do
        base=`basename $file .pdb`
        pdbsecstr $file ${base}.sec
        onehelixres.py ${base}.sec ${base}.1hr
        #rm ${base}.sec
        #find . -size 0 -delete
done

for file in *.1hr
do
        base=`basename $file .1hr`
        pdbatomselect ${base}.pdb ${base}.sel
        pdbgetresidues ${base}.1hr ${base}.sel ${base}.1res
        #rm ${base}.sel 
        #find . -size 0 -delete
done 
for file in *.1res
do
        base=`basename $file .1res`
        ca_res_organiser.py ${base}.1hr ${base}.1res ${base}.1format

for file in *.1format
do
        base=`basename $file .1format`
        #non_proline_bend_angle.py ${base}.1format ${base}.1angle
        non_proline_middle_angle.py ${base}.1format ${base}.1angle
        #rm ${base}.1format
done

"""




if __name__=="__main__":
    subprocess.call(['pdbsecstr', '../existing_code/1ct5.pdb', '../existing_code/1ct5.sec'])
    #x = subprocess.call(['pdbsecstr', '../existing_code/1ct5.pdb'])

    input_file,helix_type = core.optional_linux_argument()
    print(helix_type)
    if helix_type == 1:
        #Do the non proline 1 helix option
    if helix_type == 2:
        # Do the proline 2 helix option
    print(input_file)
    subprocess.call(['python', 'onehelixres_refactored.py', '../existing_code/1ct5.sec', '../existing_code/1ct5.1hr'])
    print('Start Main twohelixres')
    subprocess.call(['python', 'twohelixres_refactored.py', '../existing_code/1h3l.sec', '../existing_code/1h3l.2hr'])
    print('Finish Main')