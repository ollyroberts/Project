#!/bin/sh
for file in *.pdb
do
	base=`basename $file .pdb`
	pdbsecstr $file ${base}.sec
	twohelixres.py ${base}.sec ${base}.2hr
	rm ${base}.sec
	find . -size 0 -delete
done

for file in *.2hr
do
	base=`basename $file .2hr`
	pdbatomsel ${base}.pdb ${base}.sel
	pdbgetresidues ${base}.2hr ${base}.sel ${base}.res
	ca_res_organiser.py ${base}.2hr ${base}.res ${base}.format
	proline_bend_angle.py ${base}.format ${base}.angle
	rm ${base}.format
	rm ${base}.sel
	rm ${base}.2hr
	rm ${base}.res

done | less > /home/oliver/Documents/protein_counts/1hr.txt
