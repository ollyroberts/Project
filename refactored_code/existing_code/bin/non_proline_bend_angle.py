#!/usr/bin/python3
import os
import sys
import re
import numpy as np


import argparse

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
            
	if sys.platform =='linux':
		file = linux_arguments()


def windows_arguments():

	test_files = [
	('5dvi.format','5dvi.agl'),('5dz2.format','5dz2.agl'),('5e2x.format','5e2x.agl'),('5klo.format','5klo.agl'),
	('2q8g.format','2q8g.agl'),('3ix3.format','3ix3.agl'),('4beu.format','4beu.agl'),('4pmo.format','4pmo.agl'),
	('4xr8.format','4xr8.agl'),('1f0x.format','1f0x.agl'),('1mpx.format','1mpx.agl')
	]

	master('1m1j.format','1m1j.agl')
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
	

	pattern = re.compile(r'(\w+\s*?\d+?)\n(.+?)\n')
	match = pattern.findall(file_txt)
	for x in match:
		helix_name = (x[0])
		pdb_ca = x[1]



		positions = first_mid_last_finder(x[1])

		atom_inf = information_extractor(positions,pdb_ca)

		# This takes the information extracted by the information extractor and creates an atom object
		first_atom = caAtom(atom_inf[0][0],atom_inf[0][1],atom_inf[0][2],atom_inf[0][3],atom_inf[0][4])
		mid_atom = caAtom(atom_inf[1][0],atom_inf[1][1],atom_inf[1][2],atom_inf[1][3],atom_inf[1][4])
		last_atom = caAtom(atom_inf[2][0],atom_inf[2][1],atom_inf[2][2],atom_inf[2][3],atom_inf[2][4])
		
		temp = helix([first_atom,mid_atom,last_atom])
		helix_list.append(temp)

	return(helix_list)


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

	cord_p = re.compile(r'ATOM\s+?\d+?\s+?CA\s+?(\w+?)\s(\w\s*?\d+?)\s+?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)\s*?(-*?\d+?\.\d+)')

	cords 			= cord_p.findall(pdb_txt)


	first_tupple 	= (cords[first])
	middle_tupple 	= (cords[middle])
	last_tupple 	= (cords[last])

	return [first_tupple,middle_tupple,last_tupple]




def fileread(filename):
	""" opens a file object, removes the contents and strips lhwhitespace
	"""

	file = open(filename, 'r')
	output = file.read()
	output = output.lstrip()
	file.close()

	return(output)


class caAtom:
	"""  """
	atom = 'alpha carbon'

	def __init__(self, aminoacid, residue_number, x_cords , y_cords , z_cords):

		self.aa = aminoacid
		self.res_no = residue_number
		self.x = float(x_cords)
		self.y = float(y_cords)
		self.z = float(z_cords)


class helix:

	def __init__(self, atoms = []):
		self.atoms = atoms
		self.start = atoms[0].res_no

	def calculate_bond_angle(self):# self, fist_atom, second_atom, third_atom) :

		a = np.array([self.atoms[0].x,self.atoms[0].y,self.atoms[0].z])
		b = np.array([self.atoms[1].x,self.atoms[1].y,self.atoms[1].z])
		c = np.array([self.atoms[2].x,self.atoms[2].y,self.atoms[2].z])

		ba = a - b
		bc = c - b

		cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
		angle = np.arccos(cosine_angle)

		angle = (np.degrees(angle))
		self.angle = angle



def master(input_file,output_file):

	"""the purpose of this function is to create a list of residues names A11
	of proteins that are made of two helixes seperated by a gap 
	pro_eitherside is how many side of the gap I should search """
	tempstring = ""
	pdbname = str(input_file)[:4]
	helicies = helix_creator(input_file)
	for helix in (helicies):
		helix.calculate_bond_angle()
		tempstring += pdbname + " "
		tempstring += (str(helix.start) + " " + str(helix.angle) )
		tempstring += "\n"

	file = open(output_file,'w')
	file.write(tempstring)
	file.close()

	





if __name__== "__main__":
	main()


