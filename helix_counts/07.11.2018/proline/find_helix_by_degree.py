#!/usr/bin/python
import os
import sys
import re
import numpy as np


import argparse
import subprocess

"""

Changle log 
created 04.09.2018
06.09.2018 modifed to continue running if no proline is found 



The purpose of this program is take the output of all helices containing 
pdbfilename, first res, proline/mid res and angl e.g. "1h3l A170 A177 156.444"
and to list the pdblines of those angle 



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
    """ contains a battery of test files for when in a windows enviroment
    """

    #master('1h3l.format','1h3l.agl')
        


def linux_arguments():
    """
    This determins the arguments for the program when in a linux enviroment
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('angle_file', help = 'ca atoms of proteins seperated by chains')
    parser.add_argument('output_file', help = 'output file which will have the pdb name and 1st residue with ABC angle')

    #This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    angle_file = vars(args)['angle_file']
    output_file = vars(args)['output_file']

    master(angle_file,output_file)
    """
    master("proline_output.txt","specific_angles.txt")
def fileread(filename):
    """ opens a file object, removes the contents and strips lhwhitespace
    """

    file = open(filename, 'r')
    output = file.read()
    output = output.lstrip()
    file.close()

    return(output)

def master(input_file,output_file):

    """the purpose of this function is to create a list of residues names A11
    of proteins that are made of two helixes seperated by a gap 
    pro_eitherside is how many side of the gap I should search """

    helix_angles = []

    angles_file = fileread(input_file)
    angles_file =angles_file.splitlines()

    for line in angles_file:
        line_split = line.split()
        print(line_split)

        angle = line_split[3]
        angle = angle.split(".")
        angle = angle[0]

        if angle == "164":
            helix_angles.append(line)
    
    tempstring = ""
    for x in helix_angles:
        tempstring += x +"\n"
    file = open(output_file,'w')
    file.write(tempstring)
    file.close()

    #print(helix_angles)
if __name__== "__main__":
    main()
