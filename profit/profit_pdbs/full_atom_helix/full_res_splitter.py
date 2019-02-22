#!/usr/bin/python
import os
import sys
import re

import argparse


''' this is the full_res_splitter.py 
it will take as input a example.full_res file and output multiple output.txt file
 its input will be a example.full_res file 5klo.full_res containing
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

'''
this must be split into multiple files based on the how many different chains of atoms. each helix is a continuous chain of residues
everytime a residue is not 
    a) the same chain but not the next in the sequence of residues e.g. A20,A21,A75
    b) a different chain letter (e.g. (chain)A, (chain)B 
'''

def win_or_linux():



    if sys.platform =='win32':
        sys_args = windows_arguments()


    if sys.platform =='linux2':
        sys_args = linux_arguments()

    return sys_args



def linux_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('full_format_file', help = 'contains all the pdb atom details for helicies')


    #This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    input_file = vars(args)['full_format_file']

    arguments =(input_file)
    return arguments

def windows_arguments():
#    secstr_output = files.read()
    arguments =('5klo.full_res')
    return arguments


def fileread(filename):
    file = open(filename, 'r')
    output = file.read()
    file.close()
    return(output)

def filewrite(filename,input_string):

    output = open(filename,'w')
    output.write(input_string)
    output.close()

class chain: 

   def __init__(self, name):
        self.name = name
        self.residues = []  

    

if __name__== "__main__":
        arguments = win_or_linux()
        file_contents = fileread(arguments)
        filewrite('test_file.txt',file_contents)

