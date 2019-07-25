#!/usr/bin/python3
import os
import sys
import re


import argparse
import time
'''
Created by Oliver Roberts 
V1.0 - 22.02.2019
v3.0 - 23.02.2019
'''

''' this is the full_res_splitter.py 
it will take as input a example.full_res file and output multiple output.txt file
 its input will be a example.full_res file 5klo.full_res containing

'''

example_string = '''ATOM    768  O   ILE A 118     -39.704  11.044 181.540  1.00 24.59           O  
                    ATOM    769  CB  ILE A 118     -37.353  12.458 179.701  1.00 27.66           C  
                    ATOM    770  CG1 ILE A 118     -36.898  12.765 178.263  1.00 27.17           C  
                    ATOM    771  CG2 ILE A 118     -36.805  11.111 180.166  1.00 26.54           C  
                    ATOM    772  CD1 ILE A 118     -37.590  11.884 177.214  1.00 28.63           C  
                    ATOM   4110  N   SER B  76     -45.846  18.089 199.369  1.00 31.31           N  
                    ATOM   4111  CA  SER B  76     -46.902  17.558 198.515  1.00 31.76           C  
                    ATOM   4112  C   SER B  76     -47.165  16.114 198.910  1.00 36.90           C  
                    ATOM   4113  O   SER B  76     -46.334  15.478 199.565  1.00 31.79           O'''

'''
this must be split into multiple files based on the how many different chains of atoms. each helix is a continuous chain of residues
everytime a residue is not 
    a) the same chain but not the next in the sequence of residues e.g. A20,A21,A75
    b) a different chain letter (e.g. (chain)A, (chain)B 
'''

def win_or_linux():
    '''
    Determines weather the operating system is windows or linux
    '''


    if sys.platform =='win32':
        sys_args = windows_arguments()
    if sys.platform =='linux':
        sys_args = linux_arguments()
    print(sys.platform)
    return sys_args



def linux_arguments():
    '''
    Takes the arguments from the linux command line. In this case the name of the input file
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('full_format_file', help = 'contains all the pdb atom details for helicies')


    #This creates a namespace object which allows you to treat files as if they are open
    args = parser.parse_args()

    input_file = vars(args)['full_format_file']

    arguments =(input_file)
    return arguments

def windows_arguments():
    '''
    If running on  the windows version assumes that the file name will be inputed 
    '''
    arguments =('5klo.full_res')
    return arguments


def fileread(filename):
    '''
    Creates file object and closes file 
    '''
    file = open(filename, 'r')
    output = file.read()
    output = output.strip()
    file.close()
    return(output)

def filewrite(filename,input_string):
    '''
    writes a string to a target file 
    '''
    output = open(filename,'w')
    output.write(input_string)
    output.close()



class Chain(object): 
    '''
    
    '''
    chains_dict = {}

    def __init__(self, name):
        self.cname = name
        self.residues = [] 
        Chain.chains_dict[self.cname] = self



class Residue(Chain): 

    ''' Creates a residue dictionary to check if the residue is already in the dict'''
    residues_dict = {}

    def __init__(self,res_name, residue_type, residue_number,chain):
        
        self.rname = res_name
        self.rtype = residue_type
        self.rnumber = residue_number
        Residue.residues_dict[self.rname] = self

    def add_residue (self,res_name, residue_type, residue_number):

        self.residues += Residue(res_name, residue_type, residue_number)

    def add_atom(self, atom):
        self.atoms.append(atom)

class Atom(Residue): 

    atoms_dict = {}

    def __init__(self, pdb_atom_no,atom_type,x_cord,y_cord,z_cord,full_pdb_line,res_name,residue_type,residue_number,chain):
        self.anumber = pdb_atom_no
        self.type = atom_type
        self.x = float(x_cord)
        self.y = float(y_cord)
        self.z = float(z_cord)
        self.pdb_line = full_pdb_line
        super().__init__(self,res_name,residue_type,residue_number,chain)
        Atom.atoms_dict[self.anumber] = self


def chain_maker(values_dict):
        '''
        inputs: a dictionary with keys atom_number,atom_type,residue_type,chain,res_number,x_coordinates,y_coordinates,z_coordinates,pdb_line

        '''

        # atom_number = values_dict["atom_number"]
        # atom_type = values_dict["atom_type"]
        # residue_type = values_dict["residue_type"]
        chain = values_dict["chain"]
        res_number = values_dict["res_number"]
        # x_coordinates = values_dict["x_coordinates"]
        # y_coordinates = values_dict["y_coordinates"]
        # z_coordinates = values_dict["z_coordinates"]
        # pdb_line = values_dict["pdb_line"]

        res_name = "%s%s" % (chain,res_number)
        

        # if "%s" % (chain) not in Chain.chains_dict:
        #     print(" %s not in chains dict" %(chain))
        #     current_chain = Chain(chain)

        # if "%s" % (res_name) not in Residue.residues_dict:
        #     print(" %s not in residues dict" %(res_name))



def pdbline_parser(line):

    pattern = re.compile(r'ATOM\s*?(\d+)\s*?(\S+?)\s*?(\w+?)\s(\w+?)\s*?(\d+?)\s+?([-+]?\d*\.\d+|\d+)\s*?([-+]?\d*\.\d+|\d+)\s*?([-+]?\d*\.\d+|\d+)\s*?')
    values_dict = {}

    if pattern.search(line):

        hit = pattern.search(line)

        values_dict["atom_number"] = hit.group(1)
        values_dict["atom_type"] = hit.group(2)
        values_dict["residue_type"] = hit.group(3)
        values_dict["chain"] = hit.group(4)
        values_dict["res_number"] = hit.group(5)
        values_dict["x_coordinates"] = float(hit.group(6))
        values_dict["y_coordinates"] = float(hit.group(7))
        values_dict["z_coordinates"] = float(hit.group(8))
        values_dict["pdb_line"] = line
        values_dict["res_name"] = "%s%s" %(values_dict["chain"],values_dict["res_number"])

        '''creates an Atom object to store the data '''
        return(values_dict)
        

def pdb_traverser(split_pdb):


   # split_pdb = pdb_string.splitlines()
    line_counter = 0 
    current_line = pdbline_parser(split_pdb[0])
    current_chain,current_res_number = current_line["chain"],current_line["res_number"]
    breaks = []

    for line in split_pdb:

        temp_dict = pdbline_parser(line)

        previous_chain,previous_res_number = current_chain,current_res_number
        current_chain,current_res_number = temp_dict["chain"],temp_dict["res_number"]

        #print(current_res_number,previous_res_number)
        break_found = chain_break_checker(previous_chain =previous_chain, previous_res_number =previous_res_number, current_chain =current_chain, current_res_number =current_res_number)
        
        if break_found == True:
            breaks.append(line_counter) 

        #print(temp_dict["res_name"])

        line_counter += 1 


    return (breaks)
    #print(pdb_string)

def chain_break_checker(previous_chain, previous_res_number,current_chain,current_res_number):
    # for n-line in block of txt, check if 
    #((n)line chain == (n-1)line chain) AND (residue number = (residue number or (residue number-1))
    # if its not the same and a /n and then split by /n. create a list of lists 


    previous_res_number = int(previous_res_number)
    current_res_number = int(current_res_number)
    if (((current_res_number-1) == previous_res_number) or(current_res_number == previous_res_number)) and (previous_chain == current_chain):
        return False
        
    else:
        
        print("Break found. previous residue number %s chain %s and new helix current residue number %s chain %s" %(previous_res_number,previous_chain, current_res_number, current_chain))
        return True
    #print(current_res_number,previous_res_number)

    # #if (previous_chain == current_chain) and ((previous_res_number == current_chain - 1) or (previous_res_number == current_res_number)) :
        #print("same chain")


def segment_maker(break_points):
    ''' what i get 367,699,1066
        what i want [None:367],[367:699],[699:1066],[1066:None]

    '''
    first_break = None
    last_break = None

    break_points.append(last_break)

    paired_slices = []

    first_time = True 

    for position in break_points:
        if first_time == True :
            first_of_pair = None
            first_time = False
        else:
            first_of_pair = last_of_pair

        last_of_pair = position
        pair = [first_of_pair,last_of_pair]
        paired_slices.append(pair)

        first_of_pair =last_of_pair

    return(paired_slices)

def split_by_break(split_pdbs,helix_boundries):


    text_blocks = []
    for segment in helix_boundries:

        pass
        text_block = split_pdbs[segment[0]:segment[1]]
        text_blocks.append(text_block)

        #print(text_block)
    return(text_blocks)

def first_last_res_finder(text_block):

    first_line  = text_block[0]
    last_line   = text_block[-1]

    first_line_info    = pdbline_parser(first_line)
    last_line_linfo    = pdbline_parser(last_line)

    first_residue_number    = first_line_info["res_name"]
    last_residue_number     = last_line_linfo["res_name"] 
    return(first_residue_number,last_residue_number)




if __name__== "__main__":
    start = time.clock()
    arguments = win_or_linux()

    #arguments = '5klo.full_atom'

    pdb_string = fileread(arguments)
    #pdb_string = fileread('5klo.full_atom')

    pdb_string = pdb_string.splitlines()
    break_points = pdb_traverser(pdb_string)

    helix_boundries = segment_maker(break_points)
    text_blocks = split_by_break(pdb_string,helix_boundries)


    input_file_name = arguments.split('.')
    pdb_name = input_file_name[0]
    number_of_helixces = 0

    
    for block in text_blocks:
        first_res, last_res =first_last_res_finder(block)

        new_file_name = pdb_name + "_" + first_res + "_" + last_res + '.txt'
        new_file_name =  new_file_name.lower()
        
        block_string = '\n'.join(block)


        filewrite(filename=new_file_name, input_string=block_string)
        number_of_helixces += 1

    print("number of helicies from %s file is %s" %(pdb_name,number_of_helixces))
    end = time.clock()
    print("time to run fullressplitter is %f " % end)
        #filename = 
        #string = 


