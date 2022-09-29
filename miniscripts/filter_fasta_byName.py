#!/usr/bin/env python

#Script that modify transcripts fasta file from GENCODE
from Bio import SeqIO
import argparse

#Importing arguments
parser = argparse.ArgumentParser()

parser.add_argument("input", help="Input fasta files")
parser.add_argument("output", help="Output fasta files")
parser.add_argument("idFile", help="File with transcripts ID to filter")

args = parser.parse_args()

#Collect transcripts ID to pull
with open(args.idFile) as id_handle:
    wanted = set(line.rstrip("\n").split(None,1)[0] for line in id_handle)

print("... Found %i unique identifiers in %s" % (len(wanted), args.idFile))

#Pull selected transcripts
with open (args.output, "w") as output:
    for seq_record in SeqIO.parse(args.input, "fasta"):
        new_header = seq_record.name.split("|")[0].split(".")[0]

        if new_header in wanted:
            seq_record.name = seq_record.id = seq_record.description = new_header
            output.write(seq_record.format("fasta"))

            
