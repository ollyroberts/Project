#!/bin/bash
for dir in dir_001 dir_002 dir_003 dir_004 dir_005 dir_006 dir_007 dir_008 dir_009 dir_010 dir_011 dir_012 dir_013 dir_014 dir_015 dir_016 dir_017
do 
	cd ${dir}
	results=${PWD##*/}
	cp ./non_proline_output.txt ./${results}.txt
	mv ./${results}.txt ~/Documents/helix_counts/12.12.2020/non_proline/${results}.txt
	cp ./proline_output.txt ./${results}.txt
	mv ./${results}.txt ~/Documents/helix_counts/12.12.2020/proline/${results}.txt
	cd ..
done
