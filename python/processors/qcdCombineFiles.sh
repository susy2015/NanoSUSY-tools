#!/bin/bash

ls Smear_tree_* > output.txt

counter=0
files=""
filename=0

while read p; do

	files="$files $p"
	counter=$[counter+1]
	if [ $counter -eq 1000 ]
	then
		haddnano.py Smear_tree_${filename}.root ${files}
		rm -f ${files}
		counter=0
		files=""
		filename=$[filename+1]
	fi

done < output.txt

ls Smear_tree_* > output.txt

files=""

while read p; do

	files="$files $p"

done < output.txt

haddnano.py Smear_tree.root ${files}
rm -f ${files}
#rm -f output.txt

