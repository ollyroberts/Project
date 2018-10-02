This is an explanation of the shell scripts within this folder and their sub .sh and .py scripts. 

created on 30.09.2018 by Oliver Roberts 
added: 	one_angle_ext.sh
	two_angle_ext.sh
	helix_angle_combined.sh
	mutmod_pro.sh
	one_angle_mod.sh
	two_angle_mod.sh

#################
one_angle_ext.sh#
#################

Scripts used
biopy:
	pdbsecstr
	pdbatomsel
	pdbgetresidues 

my own:
	onehelixres.py 
	ca_res_organiser.py
	non_proline_middle_angle.py 

when run in a directory it takes *.pdb and creates an .angle file that contains the non-proline midpoint angle of all helicies found in that pdbfile


#################
two_angle_ext.sh#
#################

Scripts used
biopy:
	pdbsecstr
	pdbatomsel
	pdbgetresidues 

my own:
	twohelixres.py 
	ca_res_organiser.py
	proline_middle_angle.py

Similar to one_angle_ext.sh but when run in a directory it takes pdb files and creates a .angle file that contains
 the proline midpoint angle of all helicies in that pdb file.

########################
helix_angle_combined.sh#
########################

scripts used:
	one_helix_ext.sh
	two_helix_ext.sh

so is reliant upon scripts within:
	pdbsecstr
	pdbatomsel
	pdbgetresidues 
	twohelixres.py 
	ca_res_organiser.py
	proline_middle_angle.py
	

When run in full_pdb_copy directory it runs in each sub directory dir_001 - dir_017 and runs one_helix_ext.sh and 
two_helix_ext.sh

##############
mutmod_pro.sh#
##############

scripts used
bioplib:
	pdbsecstr
my own :
	pro_res_replace.py 

This program uses pro_res_replace.py to call mutmodel to replace PRO with ARG within pdbs. it uses 
"pro_res_replace.py ${base}.sec ${base}.pdb ${base}.mod" to create a .mod version of every .pdb, files with no PRO
will have similar .mob and .pdb

##################
two_angle_send.sh#
##################

no scripts used

it is a simple program, when run in a large folder containing .pdbs it will create subfolders of 100 files. 

###############################
one_helix_mod | two_helix_mod #
###############################

Scripts used
biopy:

	pdbsecstr
	pdbatomsel
	pdbgetresidues 

my own:
	onehelixres.py 
	ca_res_organiser.py
	non_proline_middle_angle.py 

versions of one_helix_ext.sh and two_helix_ext.sh except it uses .mod instead of .pdb.

