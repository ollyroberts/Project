main.py is the core ProBend script which imports functions from core_functions.py. main.py is called on the command line with format 

'main.py 1ct5.pdb -p 2'

ProBend can be called upon all PDB files in a folder using  
for pdb in *pdb; do python3 ~/Git/Project/ProBend/main.py ${pdb} -p 2; done



example_input contains a PDB file, in this case 1ct5.pdb. 
example_output contains both the files that are created with a verbose setting including .res, .format and .angle files.
.angle files contain are the final output of the bend angle of prolines from 1ct5.pdb in a whitespace seperated format.
