explains the purpose of .py scripts 

created on 30.09.2018 by Oliver Roberts 

added: 	onehelixres.py
	twohelixres.py
	ca_res_organiser.py
	non_proline_middle_angle.py
	proline_middle_angle.py
	pro_res_replace.py
	distribution.py

###############
onehelixres.py#
###############

input: .sec file , .1hr file
output: .1hr file


finds helicies of atleast length 13 (h|H) and for each pdb file outputes them in a .1hr file
seperated by /n

###############
twohelixres.py#
###############

input: .sec file , .2hr file
output: .2hr file


finds helicies of atleast length (6{H|h})PRO(6{H|h}) and for each pdb file outputes them in a .1hr file
seperated by /n

#################
ca_res_organiser#
#################

input: .1hr file, .res file , .format file 

output: .format file 

takes the Ca res from .res (in .pdb format) for the residue numbers found in .1hr.  This outputes them 
grouped by helix as .res.
  
############################
non_proline_middle_angle.py#
############################

input: 	.format file, .angle file

output:	.angle file

This calcules the bend angle for each helix within theformat file. outputs one .angle file for each .format file. 
It calculates the angle from the the mid angles in the helix. The mid angle is calculated used pdbline and creating 
3 lines around the mid residue. The midpoint of each of these lines are found and then to create the angle. 

########################
proline_middle_angle.py# 
########################

input: 	.format file, .angle file

output:	.angle file

Similar to non_proline except it finds the midmost PRO line and calcules the angle around that using pdbline 
midpoints. 

################
pro_res_replace#
################

input: ${base}.sec ${base}.mod

output: ${base}.mod

example: pro_res_replace.py ${base}.sec ${base}.mod

uses mutmodel to mutate PRO residues into ARG.

#############
distribution#
#############

 Contains a series of mathmatical models to perform on .angle catconicated files 