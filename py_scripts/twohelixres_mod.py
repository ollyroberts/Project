#!/usr/bin/python
import os
import sys
import re

import argparse
import subprocess

"""
created 07/2018
update
- v1.1 16/10/2018

updates
- v 1.1 modified double_helix_parser function to call mutmod to change the middle proline to an ALA, then pdbsecstr to check wether that changed its
helical designation



This detects whether the system is windows or linux and then calls
either windows_arguments() or linux_arguments() which controll how
a file is opened 


output is (minus the ######### with two /n at the end)
#####################
A45
A46
A47
A48
A49
A50
A51
A52
A53
A54
A55
A56
A57


#####################
"""

def win_or_linux():



    if sys.platform =='win32':
        file = windows_arguments()

            
    if sys.platform =='linux2':
        file = linux_arguments()


def windows_arguments():
#    secstr_output = files.read()
    double_helix_parser('1lxi_mod.sec','1lxi_mod.2hr')


def linux_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('sec_file', help = 'The pdbsecstr output of a pdb file that contains residues and secondary structure ')
    parser.add_argument('2hr_file', help = 'output filename which will contain double alpha helicies residues with proline')


    #This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    sec_name = vars(args)['sec_file']
    twohr_name = vars(args)['2hr_file']


    double_helix_parser(sec_name,twohr_name)


def double_helix_parser(input_file, output_file, helicies_length = 6, helix_gap = 3, pro_eitherside = 3):
    """ the purpose of this function is to create a list of residues names A11
    of proteins that are made of two helixes seperated by a gap 
    pro_eitherside is how many side of the gap I should search 
    """
    res_no_l = []       # for residue names 
    res_name_l = []     # for amino acid names
    sec_str_l = []      # for sec structure prediction

    two_helix_l = []        # contains a list aminoacids (also a list)

    # Extracts the residue no, amino acid and secstr and signs to variables
    rx_seq = re.compile(r"^(\w+?)\s+?(\w+?)\s+?(\S)", re.MULTILINE)
    text = fileread(input_file)


    # assign the matched groups in the text to the res_no_l, res_name_l and sec_str_str
    for match in rx_seq.finditer(text):
        res_no, res_name, sec_str = match.groups()

        res_no_l.append(res_no)
        res_name_l.append(res_name)
        sec_str_l += sec_str


    # creates dictionaries for each with the chain as the key
    chains_sec_str_d = keychain_value_str(res_no_l, sec_str_l)
    chains_res_no_d = keychain_value_list(res_no_l, res_no_l)
    chains_res_name_d = keychain_value_list(res_no_l, res_name_l)



    # which a Pro is found a in the res_name_d[chain] its secstr in sec_str_d is replaced with a P
    # We will then search for this P later on 

    counter = 0 
    for chain in chains_res_name_d:
        #print(chains_res_name_d[chain])
        counter = 0 
        for residue in chains_res_name_d[chain]:
            #print(chains_res_name_d[chain][counter])
            if residue == 'PRO':
                chains_sec_str_d[chain] = chains_sec_str_d[chain][:counter] + 'P' + chains_sec_str_d[chain][counter + 1:]
                #print(chains_res_no_d[chain][counter])
            counter += 1 

    # only adds if a proline is found in the gap
    # contains 2 groups, the 1st group being the whole helix and group 2 being the gap
    for x in chains_sec_str_d:
        
        regex = "([h|H]{6,}(?:.?){1}(P)(?:.?){1}[h|H]{6,})"
        p = re.compile(r"" +regex +"")

        # if one is found it prints out the residues numbers of that helix
        for match in p.finditer(chains_sec_str_d[x]):
            # adjusted to check for Proline around the gap 1 before and 1 after
            two_helix_l += [chains_res_no_d[x][ (match.start(1)) : (match.end(1)) ]]
            match_groups =(match.groups())

            # finds the location of the proline for mutation using mutmod
            pro_res =  (x + str(match.start(2)))
            print(pro_res + " :" + match.group(2))


    tempstr = ""

    for protein in two_helix_l:
        for residue in protein:
            tempstr += (residue + "\n")
        tempstr +=("\n")


    output = open(output_file, 'w')
    output.write(tempstr)
    output.close()
    #print('#####################')
    #print(tempstr)
    #print('#####################')

def fileread(filename):
    file = open(filename, 'r')
    output = file.read()
    file.close()
    return(output)

# creates a dictionary where the list key providesr first nonwhitespace is used as the
# key in this case it is the chain of the residue number. The key value is made into 
# a string associated with each chain. 
def keychain_value_str(key_provider, dict_values):
    counter = 0
    new_dict ={}
    for x in key_provider:
#   this checks if a chain exists and then adds secstr to a dictoinary 
    # of that chain letter for : secondary structutre, residue number and residue name
        if x[0] in new_dict:
            new_dict[x[0]] += dict_values[counter]
            counter += 1 
        else:
            new_dict[x[0]] = dict_values[counter]
            counter += 1 
    return new_dict

def keychain_value_list(key_provider, dict_values):
    counter = 0
    new_dict ={}
    for x in key_provider:
        if x[0] in new_dict:
            new_dict[x[0]] += [dict_values[counter]]
            counter += 1 
        else:
            new_dict[x[0]] = [dict_values[counter]]
            counter += 1 
    return new_dict

ef shell_interface(residue_pos,input_file):
    """

    The purpose of this function is to create a string calling mutmod which is sent to the command line.
    The purpose of this function is to wrap some residue chain and residue numbers with a wrapper that 
    calls the pdbline in the command line

    input: list of 3 strings        e.g ['E115 E120', 'E132 E136', 'E148 E153']

    each string is the start and end of the resdidues involved in pdbline. The first ist the starting pdbline,
    the second is the middle pdbline and the third is the end pdb line

    output: tupple of 3 strings     
    e.g. ('pdbline A52 A57 1m1j.pdb 1m1j_line_s.pdb', 'pdbline A64 A68 1m1j.pdb 1m1j_line_m.pdb', 'pdbline A74 A79 1m1j.pdb 1m1j_line_e.pdb')

    each of these are a string which will be called on the command line 

    """

    filename = input_file


    start_helix             = commandline_wrapper(residue_pos[0],filename)
    start_helix_str         = create_mutmod_string(start_helix)
 


    mid_helix               = commandline_wrapper(residue_pos[1],filename)
    mid_helix_str           = create_mutmod_string(mid_helix)
 

    #rint("end pdbline")
    #rint("residue :", str(residue_pos[2]))
    end_helix               = commandline_wrapper(residue_pos[2],filename)
    end_helix_str           = create_mutmod_string(end_helix)

    
    #print("First point :", start_pdbline_midpoint, " ", "Second point :", middle_pdbline_midpoint, " ", "Third point :", end_pdbline_midpoint)

    return(start_pdbline_midpoint,middle_pdbline_midpoint,end_pdbline_midpoint)

def commandline_wrapper(res_pair,input_name):

    """
    This function creates a wrapper for mutmod in the form 
    'mutmodel -m [c]nnn[i] (new res with -m can be repeated) ALA' 
    """
    temp_str    = ""

    temp_str    += "pdbline"
    temp_str    += " "
    temp_str    += str(res_pair)
    temp_str    += " "
    temp_str    += input_name
    temp_str    += ".pdb"

    #print(temp_str)
    return(temp_str)


def create_mutmod_string(string):
    retval=subprocess.check_output(string, shell=True)
    retval=str(retval) # Convert from byte string
    return(retval)


def main ():
        win_or_linux()

        #a parser where I can adjust the no of chains, chain length and gap length
        #double_helix_parser(file_string)

if __name__== "__main__":
        main()


