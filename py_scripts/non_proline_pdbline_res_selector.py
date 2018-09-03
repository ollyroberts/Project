#!/usr/bin/python
import os
import sys
import re
import numpy as np


import argparse
import subprocess

"""
The purpose of this program is to extract the xyz coordinates of the first (a)
 middle proline(b) and last (c) ca atom of protein

 if there is more than one proline at same dist from mid point select the lower number
"""

def main ():
	win_or_linux()

def win_or_linux():
	"""This detects whether the system is windows or linux and then calls
	either windows_arguments() or linux_arguments() which controll how
 	a file is opened """

	if sys.platform =='win32':
		file = windows_arguments()
            
	if sys.platform =='linux2':
		file = linux_arguments()


def windows_arguments():

	test_files = [
	('5dvi.format','5dvi.agl'),('5dz2.format','5dz2.agl'),('5e2x.format','5e2x.agl'),('5klo.format','5klo.agl'),
	('2q8g.format','2q8g.agl'),('3ix3.format','3ix3.agl'),('4beu.format','4beu.agl'),('4pmo.format','4pmo.agl'),
	('4xr8.format','4xr8.agl'),('1f0x.format','1f0x.agl'),('1mpx.format','1mpx.agl')
	]

	master('1h3l.format','1h3l.agl')
	#for x in test_files:
	#	master(x[0],x[1])
		

def linux_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('format_file', help = 'ca atoms of proteins seperated by chains')
	parser.add_argument('angle_file', help = 'output file which will have the pdb name and 1st residue with ABC angle')

	#This creates a namespace object which allows you to treat files as if they are open
	args = parser.parse_args()

	format_name = vars(args)['format_file']
	angle_name = vars(args)['angle_file']

	master(format_name,angle_name)



def helix_creator(input_file):
	"""
	input: string 
	the file name 

	output: list of objects 
	each object being a chain containing atom objects 

	this calls the input files and outputs a list of helicies.
	it calls first_pro_last_finder() and information_extracter() 
	to determine the 3 residues and their attributes
	"""

	file_txt = fileread(input_file)

	helix_list = []
	helix_start_mid_end = []
	temp_list = []
	

	pattern = re.compile(r'(\w+\s*?\d+?)\n(.+?)\n')
	match = pattern.findall(file_txt)
	for x in match:
		temp_list = []
		helix_name = (x[0])
		pdb_ca = x[1]



		positions = first_mid_last_finder(x[1])


		atom_inf = information_extractor(positions,pdb_ca)
		for x in atom_inf:
			temp_list += [x]
			#print(x)

		helix_start_mid_end.append(temp_list)

	return(helix_start_mid_end)


def first_mid_last_finder(protein_pdb):
	""" 
	input: a string 

	output: a tupple 
	containing 3 integers

	this finds the position of the centermost proline in the residue and
	returns what its position is in the helix
	"""
	res_p = re.compile(r'ATOM\s+?\d+?\s+?\w+?\s+?(\w+?)\s')
	res = res_p.findall(protein_pdb)



	if len(res) % 2 == 0:

		#gap window ensures only the window where prolines may be found is checked
		midpoint = int(len(res) / 2)
		pro_res = (midpoint)

	if len(res) % 2 == 1:

		midpoint = int(len(res) / 2)
		pro_res = (midpoint)

	return(0,pro_res,(len(res)-1))


def information_extractor(selected_ca,pdb_txt):
	""" 
	input: Tupple 
	containing the positions of the first/middle/last residue in the 
	aminoacid sequence 

	output:a list of 3 lists, 
	each containing the information for first/mid/last res

	This extracts the the aminoacid, resno and the xyz cords for creating the atom
	objects. The 
	"""
	first 	= selected_ca[0]
	middle 	= selected_ca[1]
	last 	= selected_ca[2]

	first_six = ""
	middle_five = ""
	last_six = ""

	cord_p = re.compile(r'ATOM\s+?\d+?\s+?CA\s+?(\w+?)\s(\w\s*?\d+?)\s+?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)')

	cords 			= cord_p.findall(pdb_txt)

	first_tupple 	= (cords[first])
	middle_tupple 	= (cords[middle])
	last_tupple 	= (cords[last])

	temp_str = (str(cords[first][1]))
	temp_str = temp_str.replace(" ", "")
	first_six += temp_str
	first_six += (" ")
	temp_str =(str(cords[first +5][1]))
	temp_str = temp_str.replace(" ", "")
	first_six += temp_str


	temp_str = (str(cords[middle -2][1]))
	temp_str =temp_str.replace(" ", "")
	middle_five += temp_str
	middle_five += (" ")
	temp_str = (str(cords[middle +2][1]))
	temp_str = temp_str.replace(" ", "")
	middle_five += temp_str


	temp_str = (str(cords[last -5][1]))
	temp_str =temp_str.replace(" ", "")
	last_six += temp_str
	last_six += (" ")
	temp_str =(str(cords[last][1]))
	temp_str = temp_str.replace(" ", "")
	last_six += temp_str

	return [first_six,middle_five,last_six]


def fileread(filename):
	""" opens a file object, removes the contents and strips lhwhitespace
	"""

	file = open(filename, 'r')
	output = file.read()
	output = output.lstrip()
	file.close()

	return(output)

def shell_interface(residue_pos,input_file):
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
	print("This helix residues are :",residue_pos)
	print("pdb file :",input_file)
	filename = input_file
	output_name = filename + ".line"

	print("start pdbline")
	start_helix = commandline_wrapper(residue_pos[0],filename,output_name)
	create_pdbline_file(start_helix)
	start_pdbline_midpoint = calculate_midpoint(output_name)

	print("middle pdbline")
	mid_helix = commandline_wrapper(residue_pos[1],filename,output_name)
	create_pdbline_file(mid_helix)
	middle_pdbline_midpoint = calculate_midpoint(output_name)

	print("end pdbline")
	print("residue :", str(residue_pos[2]))
	end_helix = commandline_wrapper(residue_pos[2],filename,output_name)
	create_pdbline_file(end_helix)
	end_pdbline_midpoint = calculate_midpoint(output_name)
	
	print("First point :", start_pdbline_midpoint, " ", "Second point :", middle_pdbline_midpoint, " ", "Third point :", end_pdbline_midpoint)

	return(start_pdbline_midpoint,middle_pdbline_midpoint,end_pdbline_midpoint)

def commandline_wrapper(res_pair,input_name,output_name):
	temp_str 	= ""

	temp_str 	+= "/home/oliver/bin/pdbline"
	temp_str 	+= " "
	temp_str 	+= str(res_pair)
	temp_str 	+= " "
	temp_str 	+= input_name
	temp_str 	+= ".pdb"
	temp_str 	+= " "
	temp_str 	+= output_name

	return(temp_str)




def create_pdbline_file(string):
	p = subprocess.Popen(string, stdout=subprocess.PIPE, shell=True)
	#p = subprocess.Popen(string, stdout=subprocess.PIPE, shell=True)
	#print("File created using :",string)


def calculate_midpoint(file_name):

	linefile = fileread(file_name)

	pdbline = re.compile(r'ATOM\s+?(\d+?)\s+?X\s+?\w+?\s\w\s*?\d+?\s+?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)')
	pdbline_xyz 	= pdbline.findall(linefile)

	

	print(" .line file findall :",pdbline_xyz)
	print("##################")
	mid = int((len(pdbline_xyz) - 1)/2)

	print(pdbline_xyz[mid])
	print(mid)

	x = pdbline_xyz[mid][1]

	y = pdbline_xyz[mid][2]

	z = pdbline_xyz[mid][3]
	print('x,y,z :',x,y,z)
	return(x,y,z)

def calculate_angle(first_xyz, second_xyz, third_xyz):

	first_xyz = list(first_xyz)
	second_xyz = list(second_xyz)
	third_xyz = list(third_xyz)

	float_first = [float(i) for i in first_xyz]
	float_second = [float(i) for i in second_xyz]
	float_third = [float(i) for i in third_xyz]


	a = np.array([float_first[0],float_first[0],float_first[0]])
	#print(a)
	b = np.array([float_second[1],float_second[1],float_second[1]])
	c = np.array([float_third[2],float_third[2],float_third[2]])

	ba = a - b
	bc = c - b

	cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
	angle = np.arccos(cosine_angle)

	angle = (np.degrees(angle))
	return(angle)

def master(input_file,output_file):

	"""the purpose of this function is to create a list of residues names A11
	of proteins that are made of two helixes seperated by a gap 
	pro_eitherside is how many side of the gap I should search """

	tempstring = ""
	pdbname = str(input_file)[:4]
	pdbline_res = helix_creator(input_file)
	print("all pdbline residues :",pdbline_res)
	for x in pdbline_res:
		one, two , three  = shell_interface(x,pdbname)
		angle =calculate_angle( one,two,three)
		print(angle)
		tempstring += input_file + " " + str(angle)

	

	#for helix in (helicies):
	#	helix.calculate_bond_angle()
	#	tempstring += pdbname + " "
	#	tempstring += (str(helix.start) + " " + str(helix.angle) )
	#	tempstring += "\n"

	file = open(output_file,'w')
	file.write(tempstring)
	file.close()

	


if __name__== "__main__":
	main()