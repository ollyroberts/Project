#!/usr/bin/python3
import os
import sys
import re
import numpy as np


import argparse
import subprocess

import proline_middle_angle as proline_mid

"""
ver 1.0 created 11.09.2018 
ver 1.1 updated on 26.09.2018 modified calculate midpoint function 

1.1 update 
-modifed taking of pdbmid point from 
	- the middle LINE residue = mid residue xyz 
to 
	- taking the first and last LINE residue xyz adding them together then halfing them
	((frist x + last x )/ 2) = midpoint x 

The purpose of this program is to extract the xyz coordinates of the first (a),
 middle proline(b), and last (c) ca atom of a protein from a .format file.

 if there is more than one proline at same dist from mid point select the lower number
 using int(len-1/2)

This calculate the bend angle using the ends of helix and the middle.
this will change the way short and long helix are reported

call structure of major functions in master function

-residue_pairs_for_pdbline -------------- first_mid_last_finder
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

The program is called on the command line, then the name of the target .format
 and then the name of the output file. e.g

 non_proline_middle_anlges.py 1h3l.format 1h3l.agl

"""


def middle_angle_linux_arguments():
	"""
	This determins the arguments for the program when in a linux enviroment
	"""


	parser = argparse.ArgumentParser()
	parser.add_argument('format_file', help = 'ca atoms of proteins seperated by chains')
	parser.add_argument('angle_file', help = 'output file which will have the pdb name and 1st residue with ABC angle')
	parser.add_argument('--pdbline', action='store_true',
						help='also creates pdbfile containing pdblines used in bend angle calculation')

	#This creates a namespace object which allows you to treat files as if they are open
	args = parser.parse_args()

	format_name 	= vars(args)['format_file']
	angle_name 		= vars(args)['angle_file']


	if args.pdbline:
		print('Pdbline creation specified')

		pdbline_option = True
		print(str(pdbline_option))
	else:
		pdbline_option = False

	return(format_name,angle_name,pdbline_option)

def fileread(filename):
	""" opens a file object, removes the contents and strips lhwhitespace
	"""

	file = open(filename, 'r')
	output = file.read()
	output = output.lstrip()
	file.close()

	return(output)

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

	# captures the first residue chain no and numer as group 1 and
	# the block of pdb txt benith as group 2
	pattern = re.compile(r'(\w+\s*?\d+?)\n(.+?)\n')
	match = pattern.findall(file_txt)
	for x in match:
		temp_list = []
		helix_name = (x[0])
		first_res_no.append(helix_name)


		pdb_ca = x[1]

		positions = first_mid_last_finder(x[1])

		atom_inf = lineinformation_extractor(positions,pdb_ca)
		for x in atom_inf:
			temp_list += [x]
			#print(x)

		helix_start_mid_end.append(temp_list)

	return(helix_start_mid_end, first_res_no)


def first_mid_last_finder(protein_pdb):
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

	mid = int((len(res) -1) / 2)

	return(mid-6,mid,mid+6)


def lineinformation_extractor(selected_ca,pdb_txt):
	"""
	input: Tupple
	containing the positions of the mid -6 /PRO/mid +6 residue in the
	aminoacid sequence (for calculating the bend angle between them)

	output:a list of 3 lists,
	each containing the information for first/mid/last res

	This extracts the the aminoacid, resno and the xyz cords for creating the atom
	objects. The
	"""
	first 	= selected_ca[0]
	middle 	= selected_ca[1]
	last 	= selected_ca[2]

	first_six 	= ""
	middle_five = ""
	last_six 	= ""

	cord_p = re.compile(r'ATOM\s+?\d+?\s+?CA\s+?(\w+?)\s(\w\s*?\d+?)\s+?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)')

	cords 			= cord_p.findall(pdb_txt)



	temp_str 	= 	(str(cords[first][1]))
	temp_str 	= 	temp_str.replace(" ", "")
	first_six 	+= 	temp_str
	first_six 	+= 	(" ")
	temp_str 	=	(str(cords[first +5][1]))
	temp_str 	= 	temp_str.replace(" ", "")
	first_six 	+= 	temp_str


	temp_str 	= 	(str(cords[middle -2][1]))
	temp_str 	=	temp_str.replace(" ", "")
	middle_five += 	temp_str
	middle_five += 	(" ")
	temp_str 	= 	(str(cords[middle +2][1]))
	temp_str 	= 	temp_str.replace(" ", "")
	middle_five += 	temp_str


	temp_str 	= 	(str(cords[last -5][1]))
	temp_str 	=	temp_str.replace(" ", "")
	last_six 	+= 	temp_str
	last_six 	+= 	(" ")
	temp_str 	=	(str(cords[last][1]))
	temp_str 	= 	temp_str.replace(" ", "")
	last_six 	+= 	temp_str
	return [first_six,middle_five,last_six]




def shell_interface(residue_pos,input_file,pdbline_option):
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



	start_helix 			= commandline_wrapper(residue_pos[0],filename)
	start_helix_str 		= create_pdbline_string(start_helix)
	start_pdbline_midpoint 	= calculate_midpoint(start_helix_str)


	mid_helix 				= commandline_wrapper(residue_pos[1],filename)
	mid_helix_str 			= create_pdbline_string(mid_helix)
	middle_pdbline_midpoint = calculate_midpoint(mid_helix_str)


	end_helix 				= commandline_wrapper(residue_pos[2],filename)
	end_helix_str 			= create_pdbline_string(end_helix)
	end_pdbline_midpoint 	= calculate_midpoint(end_helix_str)


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
	return(start_pdbline_midpoint,middle_pdbline_midpoint,end_pdbline_midpoint)

def commandline_wrapper(res_pair,input_name):
	temp_str 	= ""

	temp_str 	+= "pdbline"
	temp_str 	+= " "
	temp_str 	+= str(res_pair)
	temp_str 	+= " "
	temp_str 	+= input_name
	temp_str 	+= ".pdb"

	print(temp_str)
	return(temp_str)


def create_pdbline_string(string):
	retval=subprocess.check_output(string, shell=True)
	retval=str(retval) # Convert from byte string
	return(retval)


def calculate_midpoint(pdbline_str):

	""" 
	input: a string containing the residue information from a pdb file 

	output: a tupple of three floats  representing the x,y,z midpoints of hte pdblines 

	This calculate the bend angle by measureing pdblines created at the tips and
	the middle
	"""
	# regex extractor to find the xyz cords in the .format file 
	pdbline = re.compile(r'ATOM\s+?(\d+?)\s+?X\s+?\w+?\s\w\s*?\d+?\s+?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)')
	pdbline_xyz 	= pdbline.findall(pdbline_str)

	#finds len of number of residues 
	no_res = int(len(pdbline_xyz))


	# finds the x,y,z of the first resdue  
	first_x, first_y, first_z = pdbline_xyz[no_res -1][1], pdbline_xyz[no_res -1][2], pdbline_xyz[no_res -1][3]


	# finds the x,y,z of the last residue 
	last_x, last_y, last_z = pdbline_xyz[0][1], pdbline_xyz[0][2], pdbline_xyz[0][3]


	# finds the midpoint x,yz, of the first x,y,z and the last x,y,z 
	mid_x = ((float(first_x) + float(last_x))/2)
	mid_y = ((float(first_y) + float(last_y))/2)
	mid_z = ((float(first_z) + float(last_z))/2)

	return(mid_x, mid_y,mid_z)


def calculate_angle(first_xyz, second_xyz, third_xyz):

	first_xyz = list(first_xyz)
	second_xyz = list(second_xyz)
	third_xyz = list(third_xyz)

	float_first = [float(i) for i in first_xyz]
	float_second = [float(i) for i in second_xyz]
	float_third = [float(i) for i in third_xyz]


	a = np.array([float_first[0],float_first[1],float_first[2]])
	b = np.array([float_second[0],float_second[1],float_second[2]])
	c = np.array([float_third[0],float_third[1],float_third[2]])

	ba = a - b
	bc = c - b

	cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
	angle = np.arccos(cosine_angle)

	angle = (np.degrees(angle))
	return(angle)



	


if __name__== "__main__":


	format_file, angle_file,pdbline_option = middle_angle_linux_arguments()
	"""the purpose of this function is to create a list of residues names A11
	of proteins that are made of two helixes seperated by a gap 
	pro_eitherside is how many side of the gap I should search """

	tempstring = ""
	pdbname = str(format_file)[:4]

	pdbline_res,helix_names = residue_pairs_for_pdbline(format_file)
	print("all pdbline residues :",pdbline_res)
	print(helix_names)

	res_counter = 0
	for x in pdbline_res:


		one, two , three  = shell_interface(x,pdbname,pdbline_option)
		angle =calculate_angle( one,two,three)

		#pdbline_res takes the residue for the first pdbline and splits by space, giving the first
		first_res = str(pdbline_res[res_counter][0])
		first_res = first_res.split(" ")
		tempstring += format_file[:4] + " " + str(helix_names[res_counter]) + " " + str(angle)+ "\n"
		#print(tempstring)
		res_counter += 1

	print(tempstring)
	file = open(angle_file,'w')
	file.write(tempstring)
	file.close()

	print("Hello World")
