#!/usr/bin/python
import os
import sys
import re

import argparse


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
    double_helix_parser('1mpx.sec','1mpx.2hr')


def linux_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('sec_file', help = 'The pdbsecstr output of a pdb file that contains residues and secondary structure ')
    parser.add_argument('2hr_file', help = 'output filename which will contain double alpha helicies residues with proline')


    #This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    sec_name = vars(args)['sec_file']
    twohr_name = vars(args)['2hr_file']


    double_helix_parser(sec_name,twohr_name)

# the purpose of this function is to create a list of residues names A11
# of proteins that are made of two helixes seperated by a gap 
# pro_eitherside is how many side of the gap I should search 
def double_helix_parser(input_file, output_file, helicies_length = 6, helix_gap = 3, pro_eitherside = 3):

    res_no_l = []       # for residue names 
    res_name_l = []     # for amino acid names
    sec_str_l = []      # for sec structure prediction

    two_helix_l = []        # contains one a list aminoacids (also a list)

    # Extracts the residue no, amino acid and secstr and signs to variables
    rx_seq = re.compile(r"^(\w+?)\s+?(\w+?)\s+?(\S)", re.MULTILINE)
    text = fileread(input_file)

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
        
        regex = "(H{"+ str(helicies_length) +",}(?:([^H]{1,"+ str(helix_gap) + "})H{" + str(helicies_length) + ",}))"
        p = re.compile(r"" +regex +"")

        # if one is found it prints out the residues numbers of that helix
        for match in p.finditer(chains_sec_str_d[x]):
            # adjusted to check for Proline around the gap 1 before and 1 after
            if "PRO" in chains_res_name_d[x][ (match.start(2) -pro_eitherside) : (match.end(2)+pro_eitherside) ]:
                two_helix_l += [chains_res_no_d[x][ (match.start(1)) : (match.end(1)) ]]

    tempstr = ""

    for protein in two_helix_l:
        for residue in protein:
            tempstr += (residue + "\n")
        tempstr +=("\n")


    output = open(output_file, 'w')
    output.write(tempstr)
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


