#!/usr/bin/python
import os
import sys
import re

import argparse
import subprocess


#This detects whether the system is windows or linux and then calls
# either windows_arguments() or linux_arguments() which controll how
# a file is opened 

def win_or_linux():



    if sys.platform =='win32':
        file = windows_arguments()

            
    if sys.platform =='linux2':
        file = linux_arguments()


def windows_arguments():
#    secstr_output = files.read()
    single_helix_parser('1lxi_mod.sec','1lxi_mod.1hr')
    #ca_atom_organiser('1f0x.2hr','1f0x.res', '1f0x.format')

def linux_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('sec_file', help = 'The pdbsecstr output of a pdb file that contains residues and secondary structure ')
    parser.add_argument('1hr_file', help = 'output filename which will contain single alpha helicies residues')


    #This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    sec_name = vars(args)['sec_file']
    onehr_name = vars(args)['1hr_file']


    single_helix_parser(sec_name,onehr_name)

# the purpose of this function is to create a list of residues names A11
# of proteins that are made of two helixes seperated by a gap 
# pro_eitherside is how many side of the gap I should search 
def single_helix_parser(input_file, output_file, helicies_length = 13):

    res_no_l = []       # for residue names 
    res_name_l = []     # for amino acid names
    sec_str_l = []      # for sec structure prediction

    one_helix_l = []        # contains one a list aminoacids (also a list)

    text = fileread(input_file)

    # Extracts the residue no, amino acid and secstr and signs to variables
    rx_seq = re.compile(r"^(\w+?)\s+?(\w+?)\s+?(\S)", re.MULTILINE)

    for match in rx_seq.finditer(text):
        res_no, res_name, sec_str = match.groups()

        res_no_l.append(res_no)
        res_name_l.append(res_name)
        sec_str_l += sec_str

    chains_sec_str_d = keychain_value_str(res_no_l, sec_str_l)

    chains_res_no_d = keychain_value_list(res_no_l, res_no_l)

    chains_res_name_d = keychain_value_list(res_no_l, res_name_l)


    # only adds if a proline is found in the gap
    # contains 2 groups, the 1st group being the whole helix and group 2 being the gap
    for x in chains_sec_str_d:
        
        regex = "([H|h]{"+ str(helicies_length) +",})"
        p = re.compile(r"" +regex +"")

        # if one is found it prints out the residues numbers of that helix
        for match in p.finditer(chains_sec_str_d[x]):
            # adjusted to check for Proline around the gap 1 before and 1 after

            one_helix_l += [chains_res_no_d[x][ (match.start(1)) : (match.end(1))]]


    tempstr = ""

    for protein in one_helix_l:
        for residue in protein:
            tempstr += (residue + "\n")
        tempstr +=("\n")


    output = open(output_file, 'w')
    output.write(tempstr)
    print(tempstr)
    output.close()


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

def main ():
        win_or_linux()

        #a parser where I can adjust the no of chains, chain length and gap length
        #double_helix_parser(file_string)

if __name__== "__main__":
        main()


