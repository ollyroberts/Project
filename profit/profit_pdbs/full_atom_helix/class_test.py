#!/usr/bin/python
import os
import sys
import re

import argparse


''' this is the full_res_splitter.py 
it will take as input a example.full_res file and output multiple output.txt file
 its input will be a example.full_res file 5klo.full_res containing

'''
example_string = '''
ATOM    768  O   ILE A 118     -39.704  11.044 181.540  1.00 24.59           O  
ATOM    769  CB  ILE A 118     -37.353  12.458 179.701  1.00 27.66           C  
ATOM    770  CG1 ILE A 118     -36.898  12.765 178.263  1.00 27.17           C  
ATOM    771  CG2 ILE A 118     -36.805  11.111 180.166  1.00 26.54           C  
ATOM    772  CD1 ILE A 118     -37.590  11.884 177.214  1.00 28.63           C  
ATOM   4110  N   SER B  76     -45.846  18.089 199.369  1.00 31.31           N  
ATOM   4111  CA  SER B  76     -46.902  17.558 198.515  1.00 31.76           C  
ATOM   4112  C   SER B  76     -47.165  16.114 198.910  1.00 36.90           C  
ATOM   4113  O   SER B  76     -46.334  15.478 199.565  1.00 31.79           O  
'''




class chain: 

   def __init__(self, name):
        self.name = name
        self.residues = []  

class Atom: 

   def __init__(self, pdb_atom_no,atom_type,x_cord,y_cord,z_cord,full_pdb_line):
        self.number = pdb_atom_no
        self.type = atom_type
        self.x = float(x_cord)
        self.y = float(y_cord)
        self.z = float(z_cord)
        self.pdb_line = full_pdb_line

class Residue: 

	''' Creates a residue dictionary to check if the residue is already in the dict'''
	residues_dict = {}

	def __init__(self,residue_type, chain, residue_number):
   		self.rnumber = residue_number
   		self.rname = "%s%s" % (chain,residue_number)
   		self.rchain= chain

   		Residue.residues_dict[self.rname] = self

   	def add_atom(self, atom):
   		self.atoms.append(atom)


def pdbline_parser(list_of_pdb_lines):

	pattern = re.compile(r'ATOM\s*?(\d+)\s*?(\S+?)\s*?(\w+?)\s(\w+?)\s*?(\d+?)\s+?([-+]?\d*\.\d+|\d+)\s*?([-+]?\d*\.\d+|\d+)\s*?([-+]?\d*\.\d+|\d+)\s*?')
	atom_number = []

	#This creates a list with the index locations of each item in aa_list
	# in ca_list
	for line in list_of_pdb_lines:

		if pattern.search(line):

			hit = pattern.search(line)

			atom_number = hit.group(1)
			atom_type = hit.group(2)
			residue_type = hit.group(3)
			chain = hit.group(4)
			res_number = hit.group(5)
			x_coordinates = hit.group(6)
			y_coordinates = hit.group(7)
			z_coordinates = hit.group(8)


			pdb_line = line

			#creates an Atom object to store the data 
			atom_inf = Atom(atom_number,atom_type,x_coordinates,y_coordinates,z_coordinates,pdb_line)
			
			# checks to see if a residue of that number already exists, and if not creates it
			if "%s%s" % (chain, res_number) not in Residue.residues_dict:
				print("not in dict")
				residue_inf = Residue(residue_type =residue_type, chain=chain, residue_number=residue_type)

			for res_name in Residue.residues_dict:
				print(res_name)
				#for entry in Residue.residues_dict[res_name]:
				#	print(entry)


	
split = example_string.splitlines()
pdbline_parser(split)
#print(split)

