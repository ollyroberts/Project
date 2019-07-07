#!/usr/bin/python3
import os
import sys
import re

import argparse

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

    :param input_file:
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


def filewrite_nestedlist(filename, outerlist):
    """
    input:
    filename:name of file to be written to
    outerlist: list of lists
    This writes out the elements of a list of lists in order with a blank linke between each set of list elements and each list element on a newline
    """
    tempstr = ""

    for nestedlist in outerlist:
        for element in nestedlist:
            tempstr += (element + "\n")
        tempstr += ("\n")

    with open(filename, 'w') as out_file:
        out_file.write(tempstr)


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