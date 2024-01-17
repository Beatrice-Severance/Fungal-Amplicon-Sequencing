from Bio import SeqIO
import numpy
import os
import sys

# Specify the path to your FASTA file
fasta_file = sys.argv[1]

# Read in the FASTA file
records = list(SeqIO.parse(fasta_file, "fasta"))

# Get the base name of the file (excluding the directory and extension)
file_name = os.path.splitext(os.path.basename(fasta_file))[0]

# Modify the headers with the file name
for record in records:
    record.id = file_name
    record.description = ""

# Specify the output file name
output_file = file_name + "_newheaders" + ".fasta"

# Write the modified records to the output file
SeqIO.write(records, output_file, "fasta")

#print(f"Headers in {fasta_file} have been replaced with '{output_file}' and saved to {output_file}")
