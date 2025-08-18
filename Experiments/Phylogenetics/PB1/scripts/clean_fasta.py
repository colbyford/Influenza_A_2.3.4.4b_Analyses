#!/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import sys

# Define allowed nucleotide characters (IUPAC)
valid_nt = set("atcgryswkmbdhvn")

def clean_seq(seq):
    # Uppercase and replace invalid characters with N
    return ''.join([nt if nt in valid_nt else 'N' for nt in seq.lower()])

if len(sys.argv) < 2:
    print("Usage: python clean_fasta.py <input_fasta>")
    sys.exit(1)

input_fasta = sys.argv[1]
output_fasta = sys.argv[1].replace(".fasta", "_cleaned.fasta")

cleaned_records = []
for record in SeqIO.parse(input_fasta, "fasta"):
    cleaned_seq = clean_seq(str(record.seq))
    cleaned_records.append(SeqRecord(Seq(cleaned_seq), id=record.id, description=""))

SeqIO.write(cleaned_records, output_fasta, "fasta")
