"""
@author: jcuamatzi

@Description: this python script takes a reference genome (multifasta file) and split it into individual fasta files

This program can be used to split any other multifasta file, just is required that each sequence have is header (>)

USAGE:

python3 Split.Reference.Genome.py -f <file.fasta>

"""

from Bio import SeqIO
import argparse
import os

# Define and call arguments
parser = argparse.ArgumentParser(description="Split a multi-fasta file in individual files")
parser.add_argument("-f", "--fastaFile", help = "fasta file", required=True)
#
args = parser.parse_args() 

# Split the path and the reference file
path, reference_genome = os.path.split(args.fastaFile)

# Move to the path where is the multi fasta file
os.chdir(path)

# Open the multi-fasta file and iterate through the records
with open(reference_genome, "r") as block:
    for record in SeqIO.parse(block, "fasta"):
        header = record.id
        sequence = record.seq

        # Create a new file with the header as the filename
        output_file = f"{header}.fa"
        with open(output_file, "w") as output_block:
            SeqIO.write(record, output_block, "fasta")

print("Splitting complete.")
