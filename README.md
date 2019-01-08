# ashakit
Useful, time-saving scripts and workflows for everyday use

### Quality Check  

#### Deeptools  3.1.3

###### computeGCBias  
Prerequisites:  
**Effective genome size**: computed with `effective_genome_size.sh` script or with `unique-kmers.py` from `khmer` package. The second method should be applied when only uniquely mapped reads are used for the analysis.  

```bash
unique-kmers.py -k [MEAN READS LENGTH] [genome.fa]
```
Examples (for GRCh38 downloaded from [GENCODE](https://www.gencodegenes.org/human/release_28.html):  
Mean reads length 28bp: 2414114105  
Mean reads length 50bp: 2701262066  

**Genome in .2bit format** Can be obtained from .fa using [fatotwobit](https://anaconda.org/bioconda/ucsc-fatotwobit) programm from UCSC:  
```bash
faToTwoBit in.fa out.2bit
```
**All BAM files need to be sorted and indexed before**

Basic script:  
```bash
computeGCBias -b file.bam --effectiveGenomeSize 270126066 -g genome.2bit -o output.txt -l 50 --biasPlot plot.png 
```  

#### Qualimap 2.2.1

###### Multi-sample BamQC  
Script to run quality check for multiple BAM files  

Prerequisites:
**File describing input data**  Example script for generation a file with input data from a sample sheet containing experimental condition names and IDs for each BAM file:   
bamQC_table.R

```bash
qualimap multi-bamqc -c -d bamQC_fileslist.txt -gff regions.gtf -outdir multi_bamQC -outfile multi_bamQC.pdf -r
```

### Alignment  

#### Salmon 0.12.0  

###### Quasi-mapping mode  

**Indexing**  Run once for a particular set of reference transcripts (eg. obtained from GENCODE), -k parameter should be adjusted to the mean reads length, as it represents minimum acceptabe length for a valid match

```bash
salmon index -t gencode.v28.transcripts.fa -i gencode28_salmonindex25 --type quasi -k 25
```  

**Quantification** Code for parallel (GNU parallel) processing of all .fastq.gz samples located in ./fastq directory with mean read length = 49 bp.  

```bash
ls fastq/*_trim.fastq.gz | cut -d '/' -f 2 | cut -d '.' -f 1-5 | parallel salmon quant -i /home/JAK75/Documents/reference-genome/gencode28_salmonindex25 -l SF -r fastq/'{}'.fastq.gz -o transcript_quants/'{}' -p 4 --fldMean=49 --seqBias --gcBias --validateMappings --rangeFactorizationBins=4 --numBootstraps=1000 
```

### Useful files operations

#### Samtools 1.7

###### samtools index   
Index all .bam files in the folder. Obtained from Pierre Lindenbaum ([Biostars forum](https://www.biostars.org/p/170522/))   
```bash
ls *.bam | parallel samtools index '{}'
```	


