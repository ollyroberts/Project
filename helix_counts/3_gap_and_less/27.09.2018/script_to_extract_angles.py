#!/usr/bin/python

import os
import sys

def fileread(filename):
    file = open(filename, 'r')
    output = file.read()
    file.close()
    return(output)

def comma_seperate(filedata):
	csv_data = ""
	line_data =""
	split_data = filedata.splitlines()
	for line in split_data:
		split_line = line.split()
		line_data = ",".join(split_line)
		line_data += "\n"
		csv_data += line_data
	return(csv_data)

def writefile(filename,data):
	new_filename =(filename[0:-4]+".csv")
	file = open(new_filename,'w+')
	file.write(data)
	file.close()


filenames = ("proline.txt","non_proline.txt")


if __name__ == "__main__":
	first_data = fileread(filenames[0])
	second_data = fileread(filenames[1])

	comma_seperated_first 	= comma_seperate(first_data)
	comma_seperated_second 	= comma_seperate(second_data)

	writefile(filenames[0],comma_seperated_first)
	writefile(filenames[1],comma_seperated_second)
