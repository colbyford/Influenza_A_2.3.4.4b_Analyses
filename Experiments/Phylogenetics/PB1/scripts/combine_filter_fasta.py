#!/bin/ptyhon3

import os
import sys


def read_fasta(file_path):
    with open(file_path, 'r') as f:
        seq_id = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id:
                    yield (seq_id, ''.join(seq))
                seq_id = line
                seq = []
            else:
                seq.append(line)
        if seq_id:
            yield (seq_id, ''.join(seq))

def combine_and_filter_fasta(directory, output_file):
    seen_ids = set()
    with open(output_file, 'w') as out_f:
        for filename in os.listdir(directory):
            if filename.endswith('.fasta'):
                file_path = os.path.join(directory, filename)
                for seq_id, seq in read_fasta(file_path):
                    if seq_id not in seen_ids:
                        out_f.write(f"{seq_id}\n{seq}\n")
                        seen_ids.add(seq_id)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python combine_filter_fasta.py <input_directory> <output_fasta>")
        sys.exit(1)
    input_dir = sys.argv[1]
    output_fasta = sys.argv[2]
    combine_and_filter_fasta(input_dir, output_fasta)