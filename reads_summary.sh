#!/bin/bash
#Script to compute basic statistics for sequencing reads in FASTQ format
#Currently avaliable: mean, standard deviation  
#Multiple fastq files can be submitted

USAGE="Usage: reads_summary.sh reads1.fq/reads1.fq.gz reads2.fq/reads2.fq.gz ..."

echo "Sample    Mean    SD"

for arg in "$@"  
do
	if [ "${arg##*.}" == "gz"  ]
	then
		MEAN=$(gunzip -c $arg | awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' )
		SD=$(gunzip -c $arg | awk -v mean="$MEAN" '{if(NR%4==2) {count++; dev+=length-mean} } END{print dev/(count-1)}' )
	elif [ "${arg##*.}" == "fq" ] || [ "${arg##*.}" == "fastq" ]
	then	
		MEAN=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' $arg )
                SD=$(awk -v mean="$MEAN" '{if(NR%4==2) {count++; dev+=length-mean} } END{print dev/(count-1)}' $arg )

	else
		echo "Fastq file not recognized!"  

	fi

#echo "Sample: "${arg}" "
#echo "Mean reads length: "${MEAN}" "
#echo "Standard deviation: "${SD}" "

#echo "Sample	Mean	SD"
echo ""${arg}"	"${MEAN}"	"${SD}" "

done


  

