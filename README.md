# ashakit
Useful, time-saving scripts and workflows for everyday use

### Quality Check  

#### Deeptools  

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
	
