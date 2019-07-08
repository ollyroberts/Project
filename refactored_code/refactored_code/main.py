#!/usr/bin/python3
import subprocess
import argparse
import core_functions as core

"""
#!/bin/sh
for file in *.pdb
do
        base=`basename $file .pdb`
        pdbsecstr $file ${base}.sec
        twohelixres.py ${base}.sec ${base}.2hr
        #rm ${base}.sec
        #find . -size 0 -delete
done

for file in *.2hr
do
        base=`basename $file .2hr`
        #pdbatomselect ${base}.pdb ${base}.sel
        pdbgetresidues ${base}.2hr ${base}.sel ${base}.2res
        #rm ${base}.sel 
        #find . -size 0 -delete
done 
for file in *.2res
do
        base=`basename $file .2res`
        ca_res_organiser.py ${base}.2hr ${base}.2res ${base}.2format
        #rm ${base}.2hr 
        #find . -size 0 - delete
done

for file in *.2format
do
        base=`basename $file .2format`
        #proline_bend_angle.py ${base}.2format ${base}.angle
        proline_middle_angle.py ${base}.2format ${base}.2angle
        #rm ${base}.2format
done



"""




if __name__=="__main__":

    # calls pdbsecstr
    subprocess.call(['pdbsecstr', '../existing_code/1ct5.pdb', '../existing_code/1ct5.sec'])
    #x = subprocess.call(['pdbsecstr', '../existing_code/1ct5.pdb'])

    input_file,helix_type = core.optional_linux_argument()
    print(helix_type)
    if helix_type == 1:
        #Do the non proline 1 helix option
        #subprocess.call(
        #    ['python', 'onehelixres_refactored.py', '../existing_code/1ct5.sec', '../existing_code/1ct5.1hr'])
        print('start onehelixres')
        #input_file, output_file = core.linux_arguments()
        list_of_helices, helix_secstr, helix_resno, helix_resname = core.single_helix_parser(input_file)
        print(list_of_helices)
        core.filewrite_nestedlist(output_file, list_of_helices)
        print('finish onehelixre with %d helix ' % (len(list_of_helices)))

        #print('start onehelixres')
        #
        # list_of_helices, helix_secstr, helix_resno, helix_resname = core.single_helix_parser(input_file)
        # print(list_of_helices)
        # sec_file = str(input_file[0:-3])
        # core.filewrite_nestedlist(output_file, list_of_helices)
        # print('finish onehelixre with %d helix ' % (len(list_of_helices)))
        #


    if helix_type == 2:
        # Do the proline 2 helix option
        subprocess.call(
            ['python', 'onehelixres_refactored.py', '../existing_code/1ct5.sec', '../existing_code/1ct5.1hr'])
    print(input_file)
    subprocess.call(['python', 'onehelixres_refactored.py', '../existing_code/1ct5.sec', '../existing_code/1ct5.1hr'])
    print('Start Main twohelixres')
    subprocess.call(['python', 'twohelixres_refactored.py', '../existing_code/1h3l.sec', '../existing_code/1h3l.2hr'])
    print('Finish Main')