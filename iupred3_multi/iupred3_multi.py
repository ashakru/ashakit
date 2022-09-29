#!/usr/bin/env python

##### Script to perform disorder predictions of multiple proteins from a single FASTA file
import argparse
import os
from Bio import SeqIO

# Parse arguments
parser = argparse.ArgumentParser(description='IUPred3 Multi FASTA wrapper. Example: ./iupred3_multi.py --fasta input/someProteins.fa --outputFolder output')
parser.add_argument('--fasta', type=str, required=True, help='FASTA file with amino acid sequences')
parser.add_argument('--outputFolder', type=str, required=False, help='Output folder')

args = parser.parse_args()

# Load FASTA and run IUPred3
fasta_seq = SeqIO.parse(open(args.fasta),"fasta")

i = 0
for seq in fasta_seq:
    try:
        # Print counts
        i += 1
        if i % 1000 == 0:
            print("%s sequences processed" % (i,))

        name = seq.id
        single_fasta = args.outputFolder + "/Seq_" + name + ".fa"
        output_pred = args.outputFolder + "/IUPred3_" + name + ".txt"

        SeqIO.write(seq, single_fasta, "fasta")
        os.system("./iupred3.py -a -s no " + single_fasta + " short > " + output_pred)
        os.system("rm " + single_fasta)

    except KeyboardInterrupt:
        # User interrupt the program with ctrl+c
        exit()
