#!/usr/bin/env Rscript
#Script to prepare list of bam files for Qualimat multi-sample BamQC function
#All bam files for the analysis located in ./bam
#Usage: bamQC_table.R sampleSheet.csv bam/*_unique.bam$ output.txt 

args <- commandArgs(trailingOnly=TRUE)
#Load table
table <- read.csv(args[1], stringsAsFactors=F)[,c(1,2,5)]

#List of bam files
bamfiles <- unlist(strsplit(args[2], split="/"))[1:2]
paths <- list.files(bamfiles[1], bamfiles[2])

#Match with IDs in the sample sheet
sampleID=sapply(table[,2], function(x){paste(unlist(strsplit(x, split = "[.]"))[1:2], collapse=".")})
match=match(sampleID, sapply(paths, function(x){
					paste(unlist(strsplit(x, split="[.]"))[1:2], collapse=".")}))
table$path <- paste(bamfiles[1], paths[match], sep = "/") 
table <- table[!is.na(match),]

#Save table
output=paste(getwd(), args[3], sep="/")
write.table(table[,c(1,4,3)], file=output, quote=F, sep="\t", row.names=F, col.names=F)

print(paste0("Done! Output file: ", output)) 
