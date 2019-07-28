#!/usr/bin/python3

# -*- coding: utf-8 -*-
import os
import sys
import re
import numpy as np

import argparse
import subprocess

"""

Changle log 
v 1.0 created 04.09.2018
v 1.1 06.09.2018
v 1.2  24.07.2019

v 1.1 
modifed to continue running if no proline is found 

The purpose of this program is to extract the xyz coordinates of the first (a)
 middle proline(b) and last (c) ca atom of protein

 if there is more than one proline at same dist from mid point select the lower number
 using int(len-1/2)

This calculate the bend angle using the ends of helix and the middle.
this will change the way short and long helix are reported

call structure of major functions in master function

-residue_pairs_for_pdbline -------------- linefirst_mid_last_finder
|				  			 |
|				  			 - information_extractor
|
-shell_interface ------------ commandline_wrapper
|							 |
|							 -create_pdbline_string
|							 |
|							 - calculate_midpoint					  
|
-calculate_angle

v 1.2 

- added pdbline option to the optional linux arguments. If selected a pdbline file will be created showing the lines 
used in angle calculation

- moved the contents of the main function into the "if __name__==__main__"
- made linux arguments return instead of calling master
- modified residue_pairs_for_pdbline to also return the last residues of the whole helix
- i have not yet made use of this in the main function
-added residue_extractor_from_ca_pdb
"""


def middle_angle_linux_arguments():
    """
    This determins the arguments for the program when in a linux enviroment
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('format_file', help='ca atoms of proteins seperated by chains')
    parser.add_argument('angle_file', help='output file which will have the pdb name and 1st residue with ABC angle')
    parser.add_argument('helix_type', type=int, choices=[1, 2], default='1',
                        help='specifies if it is a proline or non proline helix')
    parser.add_argument('--pdbline', action='store_true',
                        help='also creates pdbfile containing pdblines used in bend angle calculation')

    # This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    format_name = vars(args)['format_file']
    angle_name = vars(args)['angle_file']
    helix_type = vars(args)['helix_type']

    if args.pdbline:
        print('Pdbline creation specified')

        pdbline_option = True
        print(str(pdbline_option))
    else:
        pdbline_option = False

    return (format_name, angle_name, helix_type, pdbline_option)


def fileread(filename):
    """ opens a file object, removes the contents and strips lhwhitespace
    """

    file = open(filename, 'r')
    output = file.read()
    output = output.lstrip()
    file.close()

    return (output)


def residue_pairs_for_pdbline(input_file, helix_type):
    """
    input: string
    the file name

    output: list of objects


    this calls the input files and outputs a list of helicies.
    it calls first_pro_last_finder() and information_extracter()
    to determine the 3 residues and their attributes

    it calls either proline/non_proline_segment_first_mid_last_finder depending on the
    helix type
    """

    file_txt = fileread(input_file)

    helix_start_mid_end = []

    first_res_no = []
    last_res_no = []
    whole_helix_single_residues = []
    pdbline_single_residues = []

    # match an object where the data is broken up by line breaks and chain res no linebreak
    pattern = re.compile(r'(\w+\s*?\d+?)\n(.+?)(?:(\w\s*?\d+?)\s+?(?:-*?\d+?\.\d+\s*?){1,}C\s+?)\n')
    match = pattern.findall(file_txt)

    # Each x is a helix. This helix contains the ca pdb atoms of that helix
    for helix_residues in match:
        temp_list = []
        first_res_no.append(helix_residues[0])
        last_res_no.append(helix_residues[2])
        pdb_ca = helix_residues[1]

        # has the first and last residue which will be used to
        # create the pdblines. if no proline is found that helix is skipped
        try:
            positions = proline_segment_first_mid_last_finder(helix_residues[1])
        except:

            break
        # Positions are the index locations of the helix segment used for calculating
        # bend angle made

        if helix_type == 1:
            positions = non_proline_segment_first_mid_last_finder(helix_residues[1])
            print("non proline")

        elif helix_type == 2:
            positions = proline_segment_first_mid_last_finder(helix_residues[1])
            print("proline")

        if positions == None:
            break

        atom_inf = lineinformation_extractor(positions, pdb_ca)
        for helix_residues in atom_inf:
            temp_list += [helix_residues]

        helix_start_mid_end.append(temp_list)

        whole_helix_single_residues.append(residue_extractor_from_ca_pdb(pdb_ca))

        pdbline_single_residues.append(whole_helix_single_residues[0][positions[0]:positions[2] + 1])
    # print("The index of the first residues is " +str(positions[0]) +" and "+ str(positions[2]+1))
    print("halp jesus")
    return (helix_start_mid_end, first_res_no, last_res_no, whole_helix_single_residues, pdbline_single_residues)


def proline_segment_first_mid_last_finder(protein_pdb):
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


def shell_interface(residue_pos, input_file, pdbline_option):
    """
    The purpose of this function is to wrap some residue chain and residue numbers with a wrapper that
    calls the pdbline in the command line

    input: list of 3 strings 		e.g ['E115 E120', 'E132 E136', 'E148 E153']

    each string is the start and end of the resdidues involved in pdbline. The first ist the starting pdbline,
    the second is the middle pdbline and the third is the end pdb line

    output: tupple of 3 strings
    e.g. ('pdbline A52 A57 1m1j.pdb 1m1j_line_s.pdb', 'pdbline A64 A68 1m1j.pdb 1m1j_line_m.pdb', 'pdbline A74 A79 1m1j.pdb 1m1j_line_e.pdb')

    each of these are a string which will be called on the command line

    """
    print("##################################")

    filename = input_file

    start_helix = commandline_wrapper(residue_pos[0], filename)
    start_helix_str = create_pdbline_string(start_helix)
    start_pdbline_midpoint = calculate_midpoint(start_helix_str)

    mid_helix = commandline_wrapper(residue_pos[1], filename)
    mid_helix_str = create_pdbline_string(mid_helix)
    middle_pdbline_midpoint = calculate_midpoint(mid_helix_str)

    end_helix = commandline_wrapper(residue_pos[2], filename)
    end_helix_str = create_pdbline_string(end_helix)
    end_pdbline_midpoint = calculate_midpoint(end_helix_str)

    """

        # This is added as part of the optional commandline argument --pdbline, if so it writes to a new file for every
        # is based on the horrible global variable pdbline_option which is true if argument "--pdbline" is given

        # input variables;
         residue_pos (list of lists coordiantes) [['A14 A19', 'A18 A22', 'A21 A26']
         , ['A117 A122', 'A121 A125', 'A124 A129']
         , ['A156 A161', 'A160 A164', 'A163 A168'],
          ['A199 A204', 'A203 A207', 'A206 A211']]

          filename : string "1ct5"
          start_helix_str contains the atom information from pdbline for that residue pair
        """
    if pdbline_option == True:
        print('roger roger')
        temp_res_list = []
        for pair in residue_pos:
            temp_res_list += pair.split()
        temp_filename = (str(filename) + '_' + str(temp_res_list[0] + '_' + str(temp_res_list[-1])) + '_pdbline')

        # Accesses original .pdb filename
        pdbline_file_list = []
        filename += '.pdb'
        # pdbline_file_list.append(fileread(filename))

        # this is required to remove b'string' from start/mid/end_helix_str which has been converted over when it was in byte form
        start_helix_str = start_helix_str.strip('b')
        start_helix_str = start_helix_str.strip("''")
        # start_helix_str = end_helix_str.split('\n')

        mid_helix_str = mid_helix_str.strip('b')
        mid_helix_str = mid_helix_str.strip("''")
        # mid_helix_str = end_helix_str.split('\n')

        end_helix_str = end_helix_str.strip('b')
        end_helix_str = end_helix_str.strip("''")
        # end_helix_str = end_helix_str.split('\n')
        # print(end_helix_str)

        # start_helix_str =start_helix_str.decode('utf-8')

        # this is all an attempts to combine into one list and to replace '\\n' with '\n',
        # i could have just used find and replace
        pdbline_file_list.append(start_helix_str)
        pdbline_file_list.append(mid_helix_str)
        pdbline_file_list.append(end_helix_str)
        # print(pdbline_file_list)
        pdbline_file_str = ('').join(pdbline_file_list)
        pdbline_file_list = pdbline_file_str.split('\\n')

        # writes the awakward list format to a file nicly
        with open(temp_filename, 'w') as text_file:
            text_file.writelines("%s\n" % line for line in pdbline_file_list)

    # text_file.write(str(pdbline_file_str))

    # print(temp_filename)
    return (start_pdbline_midpoint, middle_pdbline_midpoint, end_pdbline_midpoint)


def commandline_wrapper(res_pair, input_name):
    temp_str = ""

    temp_str += "pdbline"
    temp_str += " "
    temp_str += str(res_pair)
    temp_str += " "
    temp_str += input_name
    temp_str += ".pdb"

    print(temp_str)
    return (temp_str)


def create_pdbline_string(string):
    retval = subprocess.check_output(string, shell=True)
    retval = str(retval)  # Convert from byte string
    return (retval)


def calculate_midpoint(pdbline_str):
    """ This calculate the bend angle by measureing pdblines created at the tips and
    the middle
    """

    pdbline = re.compile(
        r'ATOM\s+?(\d+?)\s+?X\s+?\w+?\s\w\s*?\d+?\s+?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)')
    pdbline_xyz = pdbline.findall(pdbline_str)

    # print(" .line file findall :",pdbline_xyz)
    # print("##################")
    mid = int((len(pdbline_xyz) - 1) / 2)

    # print(pdbline_xyz[mid])
    # print(mid)

    x = pdbline_xyz[mid][1]

    y = pdbline_xyz[mid][2]

    z = pdbline_xyz[mid][3]
    # print('x,y,z :',x,y,z)
    return (x, y, z)


def calculate_angle(first_xyz, second_xyz, third_xyz):
    first_xyz = list(first_xyz)
    second_xyz = list(second_xyz)
    third_xyz = list(third_xyz)

    float_first = [float(i) for i in first_xyz]
    float_second = [float(i) for i in second_xyz]
    float_third = [float(i) for i in third_xyz]

    a = np.array([float_first[0], float_first[1], float_first[2]])
    b = np.array([float_second[0], float_second[1], float_second[2]])
    c = np.array([float_third[0], float_third[1], float_third[2]])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    angle = (np.degrees(angle))
    return (angle)


def segment_middle_res(pdb_line_res_pair):
    """
    finds the proline res by taking the pdbline res pair e.g. A174 A178,
    takes the first and adds two to get A176
    """

    res_pair = pdb_line_res_pair
    res_pair = res_pair.split()

    first_res = res_pair[0]

    first_res_no = first_res[1:]
    first_res_chain = first_res[:1]

    proline_res_no = int(first_res_no) + 2
    proline_res = (first_res_chain + str(proline_res_no))

    return (proline_res)

    def residue_extractor_from_ca_pdb(helix_pdb_info):
        helix_pdb_info


def residue_extractor_from_ca_pdb(ca_string):
    single_letter_res = []
    residue_finder = re.compile(r'ATOM\s+?\d+?\s+?\w+?\s+?(\w+?)\s')
    residues = residue_finder.findall(ca_string)

    three_letter_res_d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                          'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                          'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    for residue in residues:
        single_letter_res.append(three_letter_res_d[residue])
    single_letter_res = "".join(single_letter_res)
    return (single_letter_res)


def non_proline_segment_first_mid_last_finder(protein_pdb):
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

    mid = int((len(res) - 1) / 2)

    return (mid - 6, mid, mid + 6)


if __name__ == "__main__":

    format_file, angle_file, helix_type, pdbline_option = middle_angle_linux_arguments()
    """the purpose of this function is to create a list of residues names A11
    of proteins that are made of two helixes seperated by a gap 
    pro_eitherside is how many side of the gap I should search 
    1h3l A28 A28 A58 100.71272063185602
    """
    tempstring = ""
    pdbname = str(format_file)[:4]

    pdbline_res = residue_pairs_for_pdbline(format_file, helix_type)
    print("Helix type is %s" % helix_type)

    helix_counter = 0
    for pdbline_segment in pdbline_res[0]:
        # I would insert the extracting_ca_inf_with_first_res.py scratch here
        print("x or pdbline{0] is:")
        print(pdbline_segment)
        center_pdbline_segment = segment_middle_res(pdbline_segment[1])

        one, two, three = shell_interface(pdbline_segment, pdbname, pdbline_option)
        angle = calculate_angle(one, two, three)

        # pdbline_res takes the residue for the first pdbline and splits by space, giving the first
        whole_helix_first_res = str(pdbline_res[1][helix_counter])
        whole_helix_first_res = whole_helix_first_res.replace(" ", "")
        whole_helix_last_res = str(pdbline_res[2][helix_counter])
        whole_helix_last_res = whole_helix_last_res.replace(" ", "")

        pdbline_segment_first_res = pdbline_segment[0].split()[0]
        pdbline_segment_last_res = pdbline_segment[2].split()[1]

        # print("pdbline_segment first res is " + str(pdbline_segment_first_res) + "with last res " + pdbline_segment_last_res)

        whole_helix_single_residues = pdbline_res[3][helix_counter]
        pdbline_segment_residues = pdbline_res[4][helix_counter]

        tempstring += format_file[:4] + " " + str(whole_helix_first_res) + " " + str(whole_helix_last_res) + " " \
                      + center_pdbline_segment + " " + whole_helix_single_residues + " " + " " + str(
            angle) + " " + pdbline_segment_first_res + " " + pdbline_segment_last_res + " " + pdbline_segment_residues + " " + "\n"

        helix_counter += 1

    print(tempstring)
    file = open(angle_file, 'w')
    file.write(tempstring)
    file.close()