#!/bin/sh
for file in *.pdb
do
	base=`basename $file .pdb`
	pdbsecstr $file ${base}.sec
	onehelixres.py ${base}.sec ${base}.1hr
	rm ${base}.sec
	find . -size 0 -delete
done

for file in *.1hr
do
	base=`basename $file .1hr`
	pdbatomsel ${base}.pdb ${base}.sel
	pdbgetresidues ${base}.1hr ${base}.sel ${base}.res
	ca_res_organiser.py ${base}.1hr ${base}.res ${base}.format
	rm ${base}.sel
	rm ${base}.1hr
	rm ${base}.res

done | less > /home/oliver/Documents/protein_counts/1hr.txt
