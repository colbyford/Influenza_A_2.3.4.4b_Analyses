#!/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def clean_sequence(input_fasta, output_cleaned):
    allowed = set("ATCGatcg-")  # standard IUPAC symbols + gap
    cleaned = []
    for rec in SeqIO.parse(input_fasta, "fasta"):
        cleaned_seq = ''.join([c if c in allowed else '-' for c in str(rec.seq)])
        cleaned.append(SeqRecord(Seq(cleaned_seq.upper()), id=rec.id, description=""))

    SeqIO.write(cleaned, output_cleaned, "fasta")
    print(f"Cleaned alignment written to: {output_cleaned}")

def fasta_to_tnt(fasta_file, tnt_file):
    """
    Converts a FASTA file to TNT format.
    """
    sequences = {}  # Dictionary to store sequences
    with open(fasta_file, 'r') as fasta:
        header = None
        sequence = ""
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = sequence
                header = line[1:]
                sequence = ""
            else:
                sequence += line
        if header:
            sequences[header] = sequence

    # Write TNT file
    with open(tnt_file, 'w') as tnt:
        # TNT Header (example)
        tnt.write("# TNT format\n")
        tnt.write("# Version 1.0\n")
        tnt.write("# Sequences:\n")
        for seq_id, seq_data in sequences.items():
            # TNT Sequence Header
            seq_id = seq_id.split("|")[0]
            tnt.write(f"<{seq_id}>\n")
            tnt.write(seq_data)
            tnt.write("\n")
    
    print(f"TNT converted alignment written to: {output_cleaned}")

## Variables
os.makedirs("Trees", exist_ok=True)
input_fasta = "Alignments/PB1_H5N1_2.3.4.4.b_nucleotide_cleaned_no_stop.mafft.fasta"
output_cleaned = "Trees/PB1_H5N1_2.3.4.4.b_nucleotide_cleaned_no_stop.mafft.fasta"
output_tnt = "Trees/PB1_H5N1_2.3.4.4.b_nucleotide_cleaned_no_stop.mafft.tnt"

# Clean and convert the alignment
clean_sequence(input_fasta, output_cleaned)
fasta_to_tnt(output_cleaned, output_tnt)
