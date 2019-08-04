#!/usr/bin/python3
import time
import subprocess
import argparse


bashCommand = "python2.7 ~/Downloads/KF_err_lin/Kink_Finder.py -f ~/Git/Project/refactored_code/kink_finder_tests/1mpx.pdb -o ~/Git/Project/refactored_code/kink_finder_tests/1mpx.pdb_output/ -l '145-156' -d"
#bashCommand ="ls"
#secstr_process = subprocess.run([bashCommand], check=True, stdout=subprocess.PIPE,
#                                universal_newlines=True)

secstr_process = subprocess.run([bashCommand], shell=True,stdout=subprocess.PIPE,universal_newlines=True)
secstr_output = secstr_process.stdout

print(secstr_output)


print(end - start)