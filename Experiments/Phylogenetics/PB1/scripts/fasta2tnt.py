#!/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

## Variables
os.makedirs("Trees", exist_ok=True)
input_fasta = "Alignments/PB1_H5N1_2.3.4.4.b_nucleotide_cleaned_no_stop.mafft.fasta"
output_cleaned = "Trees/PB1_H5N1_2.3.4.4.b_nucleotide_cleaned_no_stop.mafft.fasta"
output_tnt = "Trees/PB1_H5N1_2.3.4.4.b_nucleotide_cleaned_no_stop.mafft.tnt"

## Clearning the alignment
# Remove non-standard IUPAC symbols
allowed = set("ATCGatcg-")  # standard IUPAC symbols + gap
cleaned = []
for rec in SeqIO.parse(input_fasta, "fasta"):
    cleaned_seq = ''.join([c if c in allowed else '-' for c in str(rec.seq)])
    cleaned.append(SeqRecord(Seq(cleaned_seq.upper()), id=rec.id, description=""))

SeqIO.write(cleaned, output_cleaned, "fasta")
print(f"Cleaned alignment written to: {output_cleaned}")

## Convert to TNT format
records = list(SeqIO.parse(output_cleaned, "fasta"))
nchar = len(records[0].seq)
ntax = len(records)

with open(output_tnt, "w") as out:
    out.write("xread\n")
    out.write(f"{nchar} {ntax}\n")
    for rec in records:
        rec_id = rec.id.split("|")[0]
        #out.write(f"{str(rec.seq)} {rec_id}\n")
        out.write(f"{rec_id} {str(rec.seq)}\n")
    out.write(";\n")
    out.write("proc;\n")

print(f"TNT converted alignment written to: {output_cleaned}")
