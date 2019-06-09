#!/bin/sh
for file in *.pdb
do
	base=`basename $file .pdb`
	pdbsecstr $file ${base}.sec
	onehelixres.py ${base}.sec ${base}.1hr
	#rm ${base}.sec
	#find . -size 0 -delete
done

for file in *.1hr
do
	base=`basename $file .1hr`
	pdbatomselect ${base}.pdb ${base}.sel
	pdbgetresidues ${base}.1hr ${base}.sel ${base}.1res
	#rm ${base}.sel 
	#find . -size 0 -delete
done 
for file in *.1res
do
	base=`basename $file .1res`
	ca_res_organiser.py ${base}.1hr ${base}.1res ${base}.1format
	#rm ${base}.1hr 
	#rm ${base}.1res	
	#find . -size 0 - delete
done

for file in *.1format
do
	base=`basename $file .1format`
	#non_proline_bend_angle.py ${base}.1format ${base}.1angle
	non_proline_middle_angle.py ${base}.1format ${base}.1angle
	#rm ${base}.1format
done
