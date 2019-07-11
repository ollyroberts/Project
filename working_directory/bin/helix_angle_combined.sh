#!/bin/bash
for directory in dir_001 dir_002 dir_003 dir_004 dir_005 dir_006 dir_007 dir_008 dir_009 dir_010 dir_011 dir_012 dir_013 dir_014 dir_015 dir_016 dir_017
do 	
	path=`~/Documents/full_pdb_copy`
	cd ${directory}
	echo ${directory}
	one_helix_ext.sh
	touch ./one_angle_output.txt
	cat *.angle > ./non_proline_output.txt
	rm *.res *.1hr *.format *.angle 	
	two_helix_ext.sh
	touch ./two_angle_output.txt
	cat *.angle > ./proline_output.txt
	rm *.res *.2hr *.format *.angle
	cd .. 
done
