#!/usr/bin/python3
import os
import sys
import re

import argparse

"""
Created 06.2019 by Oliver Roberts 
Changelog 
v 1.0 
v 1.1 23.07.2019


v1.1 
modified double/single_helix_res_parser 
changed input from (input_file_name) to an open file type 




functions in this module and their origin 

optional_linux_argument	    one/twohelixres.py
linux_arguments  		    one/twohelixres.py
single_helix_parser		    one/twohelixres.py
double_helix_parser		    one/twohelixres.py

ca_atom_organiser		    ca_atom_organiser
write_helix_dict		    ca_atom_organiser
first_residue_pdblines	    ca_atom_organiser
aa_chains_split			    ca_atom_organiser

filewrite_nestedlist	    ca_atom_organiser
string_nestedlist		    one/twohelixres.py
keychain_value_str		    one/twohelixres.py
keychain_value_list		    one/twohelixres.py
fileread			        one/twohelixres.py




"""

def optional_linux_argument():
    """
    Takes the first arg as a file name and the second optional argument -p to specify
    non proline (1) or proline (3) containing helices. defualt val is 1
    :return:
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'pdb_file',
        help=
        'The pdbfile that we will use for pdbsecstr '
    )
    parser.add_argument("-p", "--proline", type=int, choices=[ 1, 2],
                        default=1,
                        help="1 for [H|h] helix or 2 for helix with proline HHHPHHH for example")

    args = parser.parse_args()


    if args.proline == 2:
        print("Proline containing helices file ")

        return(args.pdb_file,args.proline)
    elif args.proline == 1:
        print("non proline containing helices ".format(args))

        return(args.pdb_file,args.proline)
    else:
        print("Incorrect -p value (should be 1 or 2) selected")
        print(args)
        return(args.pdb_file,args.proline)


def linux_arguments():
    """
    uses argparse and returns the 1st argument as the imput file, and 2nd
    argument as the output file.
    Input: .sec file made using psbsecstr on a pdbfile
    output: .1hr or .2hr file dependend on the helix parser.
    returns: input_file,output_file

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'sec_file',
        help=
        'The pdbsecstr output of a pdb file that contains residues and secondary structure '
    )
    parser.add_argument(
        '1hr_file',
        help='output filename which will contain single alpha helicies residues'
    )

    # This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    sec_name = vars(args)['sec_file']
    onehr_name = vars(args)['1hr_file']
    return (sec_name, onehr_name)


def single_helix_parser(input_file, helicies_length=13):
    """

    :param input_file string:
    :param helicies_length:
    :return:one_helix_l,chains_sec_str_d,chains_res_no_d,chains_res_name_d

    examples:
            one_helix_l (helix in list format),
       [['A8', 'A9', 'A10', 'A11', 'A12', 'A13', 'A14'], ['A116', 'A117', 'A118', 'A119', 'A120']]
            chains_sec_str_d : dict
        {'A': '-----hHHHHHHHHHHHHHHHHHHHHHHHht--EEEEEeTTthHHHHHHHHHH'}
            chains_res_no_d
        {'A': ['A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10']}
            chains_res_name_d
        {'A': ['THR', 'GLY', 'ILE', 'THR', 'TYR', 'ASP', 'GLU', 'ASP', 'ARG', 'LYS', 'THR']}

    """

    res_no_l = []  # for residue names
    res_name_l = []  # for amino acid names
    sec_str_l = []  # for sec structure prediction

    one_helix_l = []  # contains one a list aminoacids (also a list)

    #text = fileread(input_file)
    text= input_file

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
        # print(x)
        regex = "([H|h]{" + str(helicies_length) + ",})"
        p = re.compile(r"" + regex + "")

        # if one is found it prints out the residues numbers of that helix
        for match in p.finditer(chains_sec_str_d[x]):
            # adjusted to check for Proline around the gap 1 before and 1 after

            one_helix_l += [
                chains_res_no_d[x][(match.start(1)):(match.end(1))]
            ]
    return (one_helix_l, chains_sec_str_d, chains_res_no_d, chains_res_name_d)



def double_helix_parser(input_file, helicies_length=6, helix_gap=3, pro_eitherside=3):
    """ the purpose of this function is to create a list of residues names A11
    of proteins that are made of two helixes seperated by a gap
    pro_eitherside is how many side of the gap I should search

   :param input_file string:
    :param helicies_length:
    :return:two_helix_l,chains_sec_str_d,chains_res_no_d,chains_res_name_d

    examples:
            two_helix_l (helix in list format),
       [['A8', 'A9', 'A10', 'A11', 'A12', 'A13', 'A14'], ['A116', 'A117', 'A118', 'A119', 'A120']]
            chains_sec_str_d : dict
        {'A': '-----hHHHHHHHHHHHHHHHHHHHHHHHht--EEEEEeTTthHHHHHHHHHH'}
            chains_res_no_d
        {'A': ['A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10']}
            chains_res_name_d
        {'A': ['THR', 'GLY', 'ILE', 'THR', 'TYR', 'ASP', 'GLU', 'ASP', 'ARG', 'LYS', 'THR']}

    """
    res_no_l = []  # for residue names
    res_name_l = []  # for amino acid names
    sec_str_l = []  # for sec structure prediction

    two_helix_l = []  # contains one a list aminoacids (also a list)

    # Extracts the residue no, amino acid and secstr and signs to variables
    rx_seq = re.compile(r"^(\w+?)\s+?(\w+?)\s+?(\S)", re.MULTILINE)
    #text = fileread(input_file)

    text = input_file

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
        # print(chains_res_name_d[chain])
        counter = 0
        for residue in chains_res_name_d[chain]:
            # print(chains_res_name_d[chain][counter])
            if residue == 'PRO':
                chains_sec_str_d[chain] = chains_sec_str_d[chain][:counter] + 'P' + chains_sec_str_d[chain][
                                                                                    counter + 1:]
            counter += 1

            # only adds if a proline is found in the gap
    # contains 2 groups, the 1st group being the whole helix and group 2 being the gap
    for x in chains_sec_str_d:

        regex = "([h|H]{6,}(?:.?){1}P(?:.?){1}[h|H]{6,})"
        p = re.compile(r"" + regex + "")

        # if one is found it prints out the residues numbers of that helix
        for match in p.finditer(chains_sec_str_d[x]):
            # adjusted to check for Proline around the gap 1 before and 1 after
            two_helix_l += [chains_res_no_d[x][(match.start(1)): (match.end(1))]]

    return(two_helix_l,chains_sec_str_d,chains_res_no_d,chains_res_name_d)

def ca_atom_organiser(aa_chains, ca_residues, output_file):
    # opens aa_chains file and splits by blank line and line
    aa_chain_residues = fileread(aa_chains)

    first_aa_in_heli = aa_chains_split(aa_chain_residues)

    # opens file of alpha carbons and splits by line
    ca_chain_string = fileread(ca_residues)
    ca_chain_list = ca_chain_string.splitlines()

    # if there is a non value for first_residue_pdblines end the function
    # It is a Dictionary of lists, each key the first residue no of a helix, each value
    # ca atoms in that helix
    pdb_output = first_residue_pdblines(first_aa_in_heli, ca_chain_list)

    return (pdb_output)


def write_helix_dict(output_file, ca_atom_of_helix_dict):
    tempstring = ""
    output = open(output_file, 'w')

    # This gives the pdb file name at the start of the document
    # pdbfilename = aa_chains[:4]
    # tempstring += pdbfilename

    for key in sorted(ca_atom_of_helix_dict):
        # tempstring +="\n"
        # tempstring +="\n" + key +"\n"
        tempstring += key + "\n"
        for value in ca_atom_of_helix_dict[key]:
            tempstring += value
        tempstring += "\n"
    output.write(tempstring)
    # print(tempstring)

    # print(tempstring)
    output.close()


def first_residue_pdblines(aa_list, pdb_ca_list):
    """
     seachres the chain list for the starting residues, when found fills a dict
     with key the first residue in a helix and the value being the pdb file lines
     of associated residues
    """

    list_of_residues = []
    location = []
    temp_dict = {}
    temp_dict2 = {}

    # captures the residue numbers with chain for pdb format
    pattern = re.compile(r'^ATOM\s+?\d+?\s+?CA\s+?\w+?\s(\w+?\s*?\d+?)\s')
    counter = 0

    # This creates a list with the index locations of each item in aa_list
    # in ca_list
    for line in pdb_ca_list:

        if pattern.search(line):
            hit = pattern.search(line)

            list_of_residues += [hit.group(1)]

    # if the file has no fist residues in a helix an empty list is returned.
    # if a res no from aa_list is in list of residues (from pdb CA list)
    # note the pdb CA line index in location
    for x in aa_list:
        if x == '':
            return
        if x == ['']:
            return
        if x in list_of_residues:
            location.append(list_of_residues.index(x))
        else:
            return
    # This reverses the AA list (Turns [A 7, A 20, A 118] to [A 118, A20, A7]
    # It captures those CA residues up to its location and then deletes them.

    location.reverse()
    counter = 0

    for x in reversed(aa_list):
        location[counter]
        temp_dict[x] = pdb_ca_list[location[counter]:]
        del pdb_ca_list[location[counter]:]
        counter += 1

    # This creates a new dict to reverses the order of the 1st dict
    # which currently goes last amino acid to first

    for key in sorted(temp_dict):
        temp_dict2[key] = temp_dict[key]
    return (temp_dict2)


def aa_chains_split(chains):
    """
    Splits the lists of chains by empty newline and then takes the
    first amino acid from each chain for naming and organising

    """

    chains = chains.lstrip()
    chains = chains.rstrip()

    myarray = chains.split("\n\n")

    temparray = []
    temparray2 = []
    temparray3 = []

    for x in myarray:
        y = x.split("\n")
        temparray += [y]

    # temparray = list(filter(None, temparray))

    for x in temparray:
        temparray2 += [x[0]]

    for x in temparray2:

        # these spaces are needed to mimic the space in pdb files between chain
        # and res number
        if (len(x[1:])) == 1:
            y = x[0] + "   " + x[1:]
        if (len(x[1:])) == 2:
            y = x[0] + "  " + x[1:]
        if (len(x[1:])) == 3:
            y = x[0] + " " + x[1:]
        if (len(x[1:])) >= 4:
            y = x[0] + "" + x[1:]
        temparray3 += [y]

    temparray3 = sorted(temparray3)

    return (temparray3)

def filewrite_nestedlist(filename, outerlist):
    """
    input:
    filename:name of file to be written to
    outerlist: list of lists
    output: a file named filename, with a /m seperated list of lists
    This writes out the elements of a list of lists in order with a blank linke between each set of list elements and
    each list element on a newline
    """
    tempstr = ""

    for nestedlist in outerlist:
        for element in nestedlist:
            tempstr += (element + "\n")
        tempstr += ("\n")

    with open(filename, 'w') as out_file:
        out_file.write(tempstr)

def string_nestedlist(outerlist):
    """
    input:
    outerlist: list of lists
    output: returns a string with a /m seperated list of lists
    This writes out the elements of a list of lists in order with a blank linke between each set of list elements and
    each list element on a newline
    """
    tempstr = ""

    for nestedlist in outerlist:
        for element in nestedlist:
            tempstr += (element + "\n")
        tempstr += ("\n")

    return(tempstr)


def keychain_value_str(key_provider, dict_values):
    """
    creates a dictionary where the list key providesr first nonwhitespace is
     used as the key in this case it is the chain of the residue number. The
      key value is made into string associated with each chain.
    :param key_provider: takes the first character from key provider
    :param dict_values: either res_no,pdbsecstr,res_type
    :return:
    """
    counter = 0
    new_dict = {}
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
    new_dict = {}
    for x in key_provider:
        if x[0] in new_dict:
            new_dict[x[0]] += [dict_values[counter]]
            counter += 1
        else:
            new_dict[x[0]] = [dict_values[counter]]
            counter += 1
    return new_dict

def fileread(filename):
    file = open(filename, 'r')
    output = file.read()
    file.close()
    return (output)

