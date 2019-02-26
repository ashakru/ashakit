#!/usr/local/bin/Rscript

#Wrapper function at metagene package 

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      The R Script that prepares customized metagene plots using metagene package
      
      Arguments:
      --bam           - character, text file with bam files to plot, created using for example ls *.bam > bam_files.txt
      --regions       - character, BED file with regions to plot, default: NA 
      --transcripts   - character, ensembl transcripts ID (without version number) used to plot 5'UTR, CDS, 3'UTR regions, default: NA
      --gtf           - character, GTF file with gene/transcriptis models, eg. downloaded from GENCODE
      --bin_count     - numeric, 
      --normalization - character, RPM or NA
      --window_size   - numeric, window size to scale given regions
      
      Example:
      ./metagene.R --bam=\"bam_files.txt\" --regions=NA --transcriptis=\"transcripts.txt\" --bin_count=3 
      --normalization=\"RPM\" --window_size=50 \n\n")
 
  q(save="no")
}
 
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
## Arg1 default
#if(is.null(args$arg1)) {
#  ## do something
#}
 
## Arg2 default
#if(is.null(args$arg2)) {
  ## do something
#}
 
## Arg3 default
#if(is.null(args$arg3)) {
#  ## do something
#}
 
## CD.