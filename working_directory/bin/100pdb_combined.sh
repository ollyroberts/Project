one_helix_ext.sh
touch ./non_proline_output.txt
cat *.1angle > ./non_proline_output.txt
#rm *.1res *.1hr *.1format *.1angle 	
two_helix_ext.sh
touch ./proline_output.txt
cat *.2angle > ./proline_output.txt
#rm *.2res *.2hr *.2format *.2angle

