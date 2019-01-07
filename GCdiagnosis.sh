#!/bin/bash
#Script to run GC content diagnostics for all bam files (indexed and sorted) in a bam folder
#Requires GNU parallel 
#Usage: GCdiagnostic.sh [effectiveGenomeSize] [genome.2bit] [mean reads length] [threads]

mkdir GCcontent

echo "$@"

#computeGCBias form deeptools software
echo "[1] Running computeGCBias ..."

ls bam/*_unique.bam | cut -d'/' -f2 | cut -d'.' -f1-5 | parallel computeGCBias -b bam/'{}'.bam  --effectiveGenomeSize ${1} -g ${2} -o GCcontent/'{}'_GCbias.txt --biasPlot GCcontent/'{}'_GCbias.png --regionSize 500 -l ${3} -p ${4}

echo "	Done!" 
    
