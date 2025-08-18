#!/bin/python3

from Bio import SeqIO
import sys

## Variables
input_fasta = sys.argv[1]
output_fasta = sys.argv[1].replace(".fasta", "_no_stop.fasta")

with open(output_fasta, "w") as out:
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = record.seq
        if len(seq) % 3 == 0 and seq[-3:] in ["taa", "tag", "tga"]:
            record.seq = seq[:-3]
        SeqIO.write(record, out, "fasta")
