#!/usr/bin/python
import os
import sys
import re
import numpy as np
import scipy.stats as stats 

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
            
	if sys.platform =='linux2':
		file = linux_arguments()


def windows_arguments():


	master('one_helix_angle.txt','two_helix_angle.txt','distribution_output.txt')

		

def linux_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('one_helix', help = 'ca atoms of proteins seperated by chains')
	parser.add_argument('two_helix', help = 'output file which will have the pdb name and 1st residue with ABC angle')

	#This creates a namespace object which allows you to treat files as if they are open
	args = parser.parse_args()

	one_helix_file = vars(args)['one_helix']
	two_helix_file = vars(args)['two_helix']

	master(one_helix_file,two_helix_file)


def fileread(filename):
	""" opens a file object, removes the contents and strips lhwhitespace
	"""

	file = open(filename, 'r')
	output = file.read()
	output = output.lstrip()
	file.close()

	return(output)

def angle_ext(input_file):

	angle = []

	helix_file = fileread(input_file)
	helix_file = helix_file.splitlines()

	counter = 0 
	for line in helix_file:
		length =(len(line.split()))
		split_line = line.split()
		angle.append((split_line[length - 1]))

	return angle

	#print(one_helix_angle)

def is_normal(input_list):
	 Need to install scipy package using pip 

	x = stats.norm.input_list(size = 100)

	return(stats.normaltest(x))

def master(one_helix_input,two_helix_input,output_file):

	one_helix_list = angle_ext(one_helix_input)

	two_helix_list = angle_ext(two_helix_input)
	



if __name__== "__main__":
	main()


