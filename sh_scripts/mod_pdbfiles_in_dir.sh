#!/bin/bash
for directory in dir_001 dir_002 dir_003 dir_004 dir_005 dir_6 dir_007 dir_008 dir_009 dir_010 dir_011 dir_012 dir_013 dir_014 dir_015 dir_016 dir_017
do 	
	path=`~/Documents/full_pdb_copy`
	cd ${directory}
	#echo ${directory}
	mutmod_pro.sh
	rm *.sec
	cd .. 
done
