#!/usr/bin/python3

import core_functions as core


if __name__ == "__main__":


    print('start twohelixres')
    input_file, output_file = core.linux_arguments()
    list_of_helices,helix_secstr,helix_resno,helix_resname = core.double_helix_parser(input_file)
    print(list_of_helices)
    core.filewrite_nestedlist(output_file,list_of_helices)
    print("finish twohelixres with %d helix" %(len(list_of_helices)))


