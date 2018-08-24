#!/bin/sh
for file in *.pdb
do
	base=`basename $file .pdb`
	pdbsecstr $file ${base}.sec
	twohelixres.py ${base}.sec ${base}.2hr
	rm ${base}.sec
	#find . -size 0 -delete
done

for file in *.2hr
do
	base=`basename $file .2hr`
	pdbatomsel ${base}.pdb ${base}.sel
	pdbgetresidues ${base}.2hr ${base}.sel ${base}.res
	rm ${base}.sel 
	#find . -size 0 -delete
done 
for file in *.res
do
	base=`basename $file .res`
	ca_res_organiser.py ${base}.1hr ${base}.res ${base}.format
	#rm ${base}.1hr 
	#rm ${base}.res	
	#find . -size 0 - delete
done

for file in *.format
do
	base=`basename $file .format`
	proline_bend_angle.py ${base}.format ${base}.angle
	#rm ${base}.format
done
