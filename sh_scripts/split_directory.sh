#!/bin/bash
maxfilesperdir=100

# loop thought all top level directories:
while IFS= read -r -d $`\0` topleveldir
do
	#enter top level subdirectory:
	cd "$topleveldir"

	declare -i filecount=0
	declare -i direcount=0

	while IFS= read -r -d $`\0` filename 
	do 
		
		if [ "$filecount" -eq 0]
		then 
			mkdir "$dircount"
		fi

		mv "$filename" "${dircount}/"
		filecount+=1

		if [ "$filecount" -ge "$maxfilesperdir" ]
		then 
			dircount+=1
			filecount=0
		
		fi
	done < <(find -type f -print0)

	cd ..

done < <(find -mindepth 1 -maxdepth 1 - type d -print0)	
