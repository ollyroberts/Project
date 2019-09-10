#!/bin/sh
for file in *.pdb
do
	base=`basename $file .pdb`
	pdbsecstr $file ${base}.sec
	twohelixres.py ${base}.sec ${base}.2hr
	#rm ${base}.sec
	#find . -size 0 -delete
done

for file in *.2hr
do
	base=`basename $file .2hr`
	pdbatomselect ${base}.pdb ${base}.sel
	pdbgetresidues ${base}.2hr ${base}.sel ${base}.2res
	#rm ${base}.2sel 
	#find . -size 0 -delete
done 
for file in *.2res
do
	base=`basename $file .2res`
	ca_res_organiser.py ${base}.2hr ${base}.2res ${base}.2format
	#rm ${base}.2hr 
	#rm ${base}.2res	
	#find . -size 0 - delete
done

for file in *.2format
do
	base=`basename $file .2format`
	#proline_bend_angle.py ${base}.2format ${base}.2angle
	proline_middle_angle.py ${base}.2format ${base}.2angle
	#rm ${base}.2format
done
