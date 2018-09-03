#!/bin/python3
import subprocess 

command = ["ls", "-l"]
output,error  = subprocess.Popen(
                    command, universal_newlines=True,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
