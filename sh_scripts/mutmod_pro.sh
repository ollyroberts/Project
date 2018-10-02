# This program is for using the pro_res_replace which calls mutmodel to replace PRO with ARG
# in pdb files
for file in *.pdb;
do	
	base=`basename $file .pdb`
	pdbsecstr ${base}.pdb ${base}.sec
done

for file in *.sec;
do	
	base=`basename $file .sec`
	pro_res_replace.py ${base}.sec ${base}.pdb ${base}.mod
done
