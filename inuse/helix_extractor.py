#!/usr/bin/python

import argparse

def main(files):
        for f in files:
                files = open(f,'r')
                secstr_output = files.read().splitlines()


                output = secstr_output
                amino_list =[]
                secstr_str = ''
                temp = []
    
                for line in output:
                        if line is not None:
                                temp = line.split()
                                amino_list.append(temp[0])
                                secstr_str += temp[2]

                                result = []
                                counter = 0

                for x in secstr_str:
                        sss_slice =  secstr_str[counter:counter+6]
                        if sss_slice == 'HHHHHH':

                                aa0 = amino_list[counter]
                                if aa0 not in result:
                                        result.append(aa0)
    
                                aa1 = amino_list[counter+1]
                                if aa1 not in result:
                                        result.append(aa1)
    
                                aa2 = amino_list[counter+2]
                                if aa2 not in result:
                                        result.append(aa2)

                                aa3 = amino_list[counter+3]
                                if aa3 not in result:
                                        result.append(aa3)

                                aa4 = amino_list[counter+4]
                                if aa4 not in result:
                                        result.append(aa4)

                                aa5 = amino_list[counter+5]
                                if aa5 not in result:
                                        result.append(aa5)
        
                        counter += 1
                #printout = ' '.join(result)
                #print(printout)
		for x in result:
			print(x)

                

if __name__== "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('file', nargs='+', help='path to the file')
	args_namespace = parser.parse_args()
	args = vars(args_namespace)['file']
	main(args)
