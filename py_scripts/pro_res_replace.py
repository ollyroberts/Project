#!/usr/bin/python
import os
import sys
import re

import argparse
import subprocess


"""
pro_res_replace.py ${base}.sec ${base}.mod

The purpose of this is to take as input : ${base}.sec ${base}.mod
example (pro_res_replace.py ${base}.sec ${base}.mod)

and then call 
"mutmodel -m A28 PRO -m A29 -m A30 PRO base.pdb base.mod"

find_prolines

pro_location
command_line_string
send_to_cmdline

"""
def win_or_linux():
    """
    detects wheather the system is windows of linux and then calls either 
    awindows_arguments() or linux_arguments() which controll how a file is opened
    """

    if sys.platform =='win32':
        file = windows_arguments()
        return file

    if sys.platform =='linux2':
        file = linux_arguments()
        return file


def windows_arguments():
#    secstr_output = files.read()
    return('1h3l.secstr','1h3l.mod')
    #ca_atom_organiser('1f0x.2hr','1f0x.res', '1f0x.format')

def linux_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('sec_file', help = 'The pdbsecstr output of a pdb file that contains residues and secondary structure used for findnig proline ')
    parser.add_argument('pdb_file', help = ' The pdbfile which will have have its residues modified')
    parser.add_argument('mod_file', help = ' The output modified file ')

    #This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    sec_name = vars(args)['sec_file']
    pdb_name = vars(args)['pdb_file']
    mod_name = vars(args)['mod_file']

    return(sec_name,pdb_name,mod_name)

def fileread(filename):
    """
    input: a string of file name 
    ooutput: a string of that object 
    """
    file = open(filename, 'r')
    output = file.read()
    file.close()
    return(output)


def find_prolines(sec_file):
    """
    input: string
    output: 2 lists, first containing a strings of res numbers, second containing res names

    This finds the locations of prolines in the .sec file (output of pdbsecstr 2)
    """

    pro_res = []
    sec_string = fileread(sec_file)
    res = sec_string.splitlines()

    for line in res:
        line_split = line.split()
        if line_split[1] == "PRO":
            pro_res.append(line_split[0])
        #line.split(" ")
        #print(line_split)

    return (pro_res)

def command_line_string(res_list, pdb_filename, output_file):
    """
    input = list of residues showing proline locations, input filename

    output = string that will be entered on the linux command line 
            eg. "mutmodel -m A28 PRO -m A29 -m A30 PRO 1h3l_copy.pdb 1h3l_copy.mod"
    """


    temp_str = ""
    temp_str += "mutmodel "

    for x in res_list:
        temp_str += "-m " + str(x) + " ARG "
    temp_str += str(pdb_filename) + " "
    temp_str += str(output_file) + " "

    if not res_list:
        temp_str ="cp " + str(pdb_filename) + " " + str(output_file)
        print("no prolines replaced")
    else:
        print("Prolines found :")
        print(res_list)

    print(temp_str)

    return(temp_str)


def send_to_cmdline(string):
    """
    uses the subprocess module to interact with linux shell and creates a new 1hr file. 

    input: string that is passed to command line 

    output: None 
    The new .mod file is created by the mutmodel on the command line 
    """
    retval=subprocess.check_output(string, shell=True)
    retval=str(retval) # Convert from byte string


def main ():
    input_output =win_or_linux()

    sec_file, pdb_file, mod_file = input_output[0],input_output[1],input_output[2]

    res_no_list = find_prolines(sec_file)
    linux_string =(command_line_string(res_no_list,pdb_file, mod_file))
    send_to_cmdline(linux_string)


if __name__== "__main__":
        main()