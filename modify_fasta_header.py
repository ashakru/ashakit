#!/usr/bin/env python

#Script that modify transcripts fasta file from GENCODE
from Bio import SeqIO
import argparse

#Importing arguments
parser = argparse.ArgumentParser()

parser.add_argument("input", help="Input fasta files")
parser.add_argument("output", help="Output fasta files")

args = parser.parse_args()

#Modify FASTA header

with open (args.output, "w") as output:
    for seq_record in SeqIO.parse(args.input, "fasta"):
        new_header = seq_record.name.split("|")[0].split(".")[0]
        seq_record.name = seq_record.id = seq_record.description = new_header
        output.write(seq_record.format("fasta")) 
