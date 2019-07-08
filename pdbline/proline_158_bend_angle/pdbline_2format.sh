for file in *.2format; do base=`basename $file .2format`;proline_middle_angle.py ${base}.2format ${base}.2angle --pdbline; done

