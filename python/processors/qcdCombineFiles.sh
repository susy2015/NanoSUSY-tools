#!/bin/bash

ls *.root > output.txt

counter=0
totcounter=0
files=""
filename=0
loopnum=0
getFinalName=''
totLines=$(wc -l output.txt | awk '{ print $1 }')

while [ $totLines -ne 1 ]
do

	loopnum=$[loopnum+1]
	
	while read p; do
	
		if [[ ${p} != *"Smear_tree"* ]]; then
			getFinalName=${p%.*}
		fi
		
		files="$files $p"
		counter=$[counter+1]
		totcounter=$[totcounter+1]
		
		if [ $counter -eq 1000 ] || [ $totcounter -eq $totLines ]
		then
			haddnano.py Smear_tree_buff_${loopnum}_${filename}.root ${files}
			rm -f ${files}
			counter=0
			files=""
			filename=$[filename+1]
		fi
	
	done < output.txt

	filename=0
	totcounter=0
	ls *.root > output.txt
	totLines=$(wc -l output.txt | awk '{ print $1 }')

done

mv Smear_tree_buff_${loopnum}_0.root ${getFinalName}_smear.root

