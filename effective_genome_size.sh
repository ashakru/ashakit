#!/bin/bash

#Script that computes effective genome size (portion of mappable genome) from provided BED or GTF file
#Prerequisites: bedtools (tested on v2.17.0) 

#GTF to BED transformation (if needed): gtf2bed function from BEDOPS package

USAGE="Usage:effective_genome_size.sh [BED file] \n -h/--help info"

for arg in "$@"
do		

	if [ "$arg" == "--help" ] || [ "$arg" == "-h" ]
	then
        	echo "$USAGE"
	elif [ "${arg##*.}" == "bed" ]
	then
		BED=${arg}
	else
		echo "$USAGE"
    fi
done

#Compute effective genome size
echo "[1] Sorting and merging ..."
SORTED=$(sort -k 1,1 -k2,2n $BED | bedtools merge)
echo "[2] Computing effective genome size ..."
SIZE=$(echo "${SORTED}" | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
echo "Effective Genome Size: $SIZE"
