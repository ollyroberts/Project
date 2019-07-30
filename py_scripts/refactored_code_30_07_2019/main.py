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
def onehelix_res(file_basename,secstr_output):

    # This if a helix with no gaps and functions like the old one_helix_res
    # akin to onehelixres.py

    filename_hr = (file_basename + ".1hr")

    filename_format = (file_basename + ".1format")

    print('start onehelixres')
    list_of_helices, helix_secstr, helix_resno, helix_resname = core.single_helix_parser(secstr_output)
    core.filewrite_nestedlist(filename_hr, list_of_helices)
    print('finish onehelixre with %d helix ' % (len(list_of_helices)))

    # subprocess.run(['pdbgetresidues', filename_hr, sel_filename, filename_res])
    # dict_of_ca_in_helix = core.ca_atom_organiser(filename_hr, filename_res, filename_format)
    # core.write_helix_dict(filename_format, dict_of_ca_in_helix)
    return(list_of_helices,file_basename,sec_filename,sel_filename)

if __name__=="__main__":

    pdb_filename = None
    sec_filename = None
    sel_filename = None

    file_basename = None
    helix_type = 1

    verbose = False


    pdb_filename, helix_type = core.optional_linux_argument()

    file_basename = pdb_filename.split('.')[0]
    sec_filename = (file_basename + ".sec")
    sel_filename = (file_basename + ".sel")


    secstr_process = subprocess.run(['pdbsecstr', pdb_filename], check=True, stdout=subprocess.PIPE, universal_newlines=True)
    secstr_output = secstr_process.stdout

    atomsel_process = subprocess.run(['pdbatomselect', pdb_filename,sel_filename])
    #atomsel_output = atomsel_process.stdout



    if helix_type == 1:
         # This if a helix with no gaps and functions like the old one_helix_res
         # akin to onehelixres.py

        filename_hr = (file_basename + ".1hr")
        filename_res = (file_basename + ".1res")
        filename_format = (file_basename + ".1format")

        print('start onehelixres')
        list_of_helices, helix_secstr, helix_resno, helix_resname = core.single_helix_parser(secstr_output)
        core.filewrite_nestedlist(filename_hr, list_of_helices)
        print('finish onehelixre with %d helix ' % (len(list_of_helices)))

        # subprocess.run(['pdbgetresidues', filename_hr, sel_filename, filename_res])
        # dict_of_ca_in_helix = core.ca_atom_organiser(filename_hr, filename_res, filename_format)
        # core.write_helix_dict(filename_format, dict_of_ca_in_helix)


    if helix_type == 2:
        # This requires the helix to have a proline
        # akin to twohelixres.py and ca_res_organiser

        filename_hr = (file_basename + ".2hr")
        filename_res = (file_basename + ".2res")
        filename_format = (file_basename + ".2format")

        print('start twohelixres')
        list_of_helices, helix_secstr, helix_resno, helix_resname = core.double_helix_parser(secstr_output)
        core.filewrite_nestedlist(filename_hr, list_of_helices)
        print("finish twohelixres with %d helix" % (len(list_of_helices)))

        # subprocess.run(['pdbgetresidues', filename_hr, sel_filename, filename_res])
        # dict_of_ca_in_helix = core.ca_atom_organiser(filename_hr, filename_res, filename_format)
        # try :
        #     core.write_helix_dict(filename_format, dict_of_ca_in_helix)
        # except:
        #     print("No proline helix found")
        #core.write_helix_dict(filename_format, dict_of_ca_in_helix)

    subprocess.run(['pdbgetresidues', filename_hr, sel_filename, filename_res])
    dict_of_ca_in_helix = core.ca_atom_organiser(filename_hr, filename_res, filename_format)
    list_of_helices_filestring = core.string_nestedlist(list_of_helices)
    print(dict_of_ca_in_helix)
    subprocess.run(['./proline_and_nonproline_middle_angle.py', filename_format, filename_format, str(helix_type)])


    """
    This segment is the ca_res_organiser_refactored
    """


    """
    pdbatomselect ${base}.pdb ${base}.sel
    dbgetresidues 1ct5.1hr 1ct5.sel 1ct5.1res
    """


    print("YAY")