#!/usr/bin/python3
import re
import numpy as np

import argparse
import subprocess

"""
Current output 
<class 'tuple'>: ([['A52 A57', 'A56 A60', 'A59 A64'], 
['B52 B57', 'B56 B60', 'B59 B64']], ['A  28', 'B  28'])
"""


def residue_pairs_for_pdbline(input_file):
    """
    input: string
    the file name

    output: list of objects


    this calls the input files and outputs a list of helicies.
    it calls first_pro_last_finder() and information_extracter()
    to determine the 3 residues and their attributes
    """

    file_txt = fileread(input_file)

    helix_start_mid_end = []

    first_res_no = []
    last_res_no = []

    # match an object where the data is broken up by line breaks and chain res no linebreak
    pattern = re.compile(r'(\w+\s*?\d+?)\n(.+?)(?:(\w\s*?\d+?)\s+?(?:-*?\d+?\.\d+\s*?){1,}C\s+?)\n')
    #pattern = re.compile(r'ATOM\s+?\d+?\s+?CA\s+?(\w+?)\s(\w\s*?\d+?)\s')
    match = pattern.findall(file_txt)

    # Each x is a helix. This helix contains the ca pdb atoms of that helix
    for x in match:
        temp_list = []
        first_res_no.append(x[0])
        last_res_no.append(x[2])
        pdb_ca = x[1]

        # has the first and last residue which will be used to
        # create the pdblines. if no proline is found that helix is skipped
        try:
            positions = linefirst_mid_last_finder(x[1])
        except:

            break

        positions = linefirst_mid_last_finder(x[1])

        if positions == None:
            break

        atom_inf = lineinformation_extractor(positions, pdb_ca)
        for x in atom_inf:
            temp_list += [x]

        helix_start_mid_end.append(temp_list)

    return (helix_start_mid_end, first_res_no, last_res_no)


def linefirst_mid_last_finder(protein_pdb):
    """
    input: a string

	output: a tupple
	containing 3 integers

	this finds the position that will be used for pdbline for creating
	lines of best fit. The first residue of the first pdbline, the proline
	at the middle of the second pdbline and the last residue of the third
	pdb line.

	"""
    res_p = re.compile(r'ATOM\s+?\d+?\s+?\w+?\s+?(\w+?)\s')

    res = res_p.findall(protein_pdb)

    counter = 0
    gap_window = int(0.5 * len(res)) - 6
    pro_res = None

    while counter <= (gap_window):

        midpoint = int((len(res) - 1) / 2)

        if res[midpoint - counter] == "PRO":

            pro_res = (midpoint - counter)
            counter += 1
            break

        elif res[midpoint + counter] == "PRO":

            pro_res = (midpoint + counter)
            counter += 1
            break

        else:

            counter += 1

    if pro_res == None:
        print("no PRO res found")
        return (None)

    return (pro_res - 6, pro_res, pro_res + 6)


def lineinformation_extractor(selected_ca, pdb_txt):
    """
    input: Tupple
    containing the positions of the mid -6 /PRO/mid +6 residue in the
    aminoacid sequence (for calculating the bend angle between them)

    output:a list of 3 lists,
    each containing the information for first/mid/last res

    This extracts the the aminoacid, resno and the xyz cords for creating the atom
    objects. The
    """
    first = selected_ca[0]
    middle = selected_ca[1]
    last = selected_ca[2]

    first_six = ""
    middle_five = ""
    last_six = ""

    cord_p = re.compile(
        r'ATOM\s+?\d+?\s+?CA\s+?(\w+?)\s(\w\s*?\d+?)\s+?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)')

    cords = cord_p.findall(pdb_txt)

    temp_str = (str(cords[first][1]))
    temp_str = temp_str.replace(" ", "")
    first_six += temp_str
    first_six += (" ")
    temp_str = (str(cords[first + 5][1]))
    temp_str = temp_str.replace(" ", "")
    first_six += temp_str

    temp_str = (str(cords[middle - 2][1]))
    temp_str = temp_str.replace(" ", "")
    middle_five += temp_str
    middle_five += (" ")
    temp_str = (str(cords[middle + 2][1]))
    temp_str = temp_str.replace(" ", "")
    middle_five += temp_str

    temp_str = (str(cords[last - 5][1]))
    temp_str = temp_str.replace(" ", "")
    last_six += temp_str
    last_six += (" ")
    temp_str = (str(cords[last][1]))
    temp_str = temp_str.replace(" ", "")
    last_six += temp_str
    return [first_six, middle_five, last_six]


def fileread(filename):
    """ opens a file object, removes the contents and strips lhwhitespace
    """

    file = open(filename, 'r')
    output = file.read()
    output = output.lstrip()
    file.close()

    return (output)


if __name__ == "__main__":
    format_file = "1h3l.2format"
    x = residue_pairs_for_pdbline(format_file)
    print(x)
