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
    ca_atom_organiser('1xzw.1hr','1xzw.res', '1xzw.format')

def linux_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('2hr_file', help = 'The output of the twohelixextractor file that that contains the residue no of the helixes is required input')
	parser.add_argument('res_file', help = 'The .res file containing the alpha carbon residues that is required input')
	parser.add_argument('formated', help = 'output file in $file.format')

	#This creates a namespace object which allows you to treat files as if they are open
	args = parser.parse_args()

	twohe_name = vars(args)['2hr_file']
	res_name = vars(args)['res_file']
	output_name = vars(args)['formated']

	ca_atom_organiser(twohe_name,res_name,output_name)

# the purpose of this function is to create a list of residues names A11
# of proteins that are made of two helixes seperated by a gap 
# pro_eitherside is how many side of the gap I should search 
def ca_atom_organiser(aa_chains,ca_residues,output_file):



	tempstring = ""
	# opens aa_chains file and splits by blank line and line 
	aa_chain_residues = fileread(aa_chains)



	first_aa_in_heli = aa_chains_split(aa_chain_residues)


	# opens file of alpha carbons and splits by line 
	ca_chain_string = fileread(ca_residues)
	ca_chain_list = ca_chain_string.splitlines()



	# if there is a non value for first_residue_pdblines end the function
	pdb_output =first_residue_pdblines(first_aa_in_heli,ca_chain_list)
	if  not pdb_output:
		return

	output = open(output_file,'w')

	# This gives the pdb file name at the start of the document 
	#pdbfilename = aa_chains[:4]
	#tempstring += pdbfilename

	for key in sorted(pdb_output):
		#tempstring +="\n"
		#tempstring +="\n" + key +"\n"
		tempstring +=key +"\n"
		for value in pdb_output[key]:
			tempstring += value
		tempstring += "\n"
	output.write(tempstring)

	#print(tempstring)
	output.close()



def first_residue_pdblines(aa_list,pdb_ca_list):
	"""
	 seachres the chain list for the starting residues, when found fills a dict
	 with key the first residue in a helix and the value being the pdb file lines 
	 of associated residues 
	"""




	list_of_residues = []
	location= []
	temp_dict = {}
	temp_dict2 = {}

	# captures the residue numbers with chain for pdb format
	pattern = re.compile(r'^ATOM\s+?\d+?\s+?CA\s+?\w+?\s(\w+?\s*?\d+?)\s')
	counter = 0

	#This creates a list with the index locations of each item in aa_list
	# in ca_list
	for line in pdb_ca_list:
		if hit == True:
			hit = pattern.search(line)
			list_of_residues += [hit.group(1)]

# if the file has no fist residues in a helix an empty list is returned.
# end program
	for x in aa_list:
		if x == '':
			return
		if x == ['']:
			return
		else:
			location.append(list_of_residues.index(x))

	location.reverse()
	counter = 0 

	for x in reversed(aa_list):
		location[counter]
		temp_dict[x] = pdb_ca_list[location[counter]:]
		del pdb_ca_list[location[counter]:]
		counter += 1 

	#This creates a new dict to reverses the order of the 1st dict
	# which currently goes last amino acid to first 

	for key in sorted(temp_dict):
		temp_dict2[key] = temp_dict[key]


	return(temp_dict2)
	


def aa_chains_split(chains):
	"""
	Splits the lists of chains by empty newline and then takes the 
	first amino acid from each chain for naming and organising

	"""

	chains = chains.lstrip()
	chains = chains.rstrip()

	myarray = chains.split("\n\n")

	temparray= []
	temparray2= []
	temparray3 = []

	for x in myarray:
		y = x.split("\n")
		temparray += [y]

	#temparray = list(filter(None, temparray))

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

	return(temparray3)




def fileread(filename):
	file = open(filename, 'r')
	output = file.read()
	file.close()
	return(output)


def main ():
	win_or_linux()



if __name__== "__main__":
	main()


