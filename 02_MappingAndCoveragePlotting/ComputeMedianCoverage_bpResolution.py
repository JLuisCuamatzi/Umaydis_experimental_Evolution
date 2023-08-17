"""
@ author: jcuamatzi

@ Description: estimate the normalized coverage at base-pair resolution

USAGE:

python3 ComputeMedianCoverage_bpResolution.py

"""

import os
import sys
import time
import pandas as pd
import gzip
import argparse

# Functions
# Function to check if the Input File Exists
def check_file_exists(file_path):
    if os.path.isfile(file_path):
        print(f"Input file for sample {sample_name} exists. The script will continue.")
    else:
        sys.exit(f"The input file {file_path} does not exist. Bye.")

# Function to Calculate the Normalized Coverage at bp-Resolution
def process_sample(input_file):
    # Get sample name from input file    
    # Open and read the input file
    with gzip.open(input_file, 'rt') as file:
        df = pd.read_csv(file, sep='\t', header=None, names=column_names)
    
    # Calculate the coverage median in all genome
    coverage_median = df['DepthCoverage'].median()
    
    # Filter by selecting chr 9
    df = df[df['Chr'] == args.chromosome]
    
    # Estimate Normalized coverage in each position
    df['normalized.coverage'] = df['DepthCoverage'] / coverage_median
    df['Sample'] = sample_name
    
    # Save the DataFrame as a text file
    # Output filename
    output_path = "normalizedCoverageTables/bpResolution/"
    os.makedirs(output_path, exist_ok=True)
    output_file = os.path.join(output_path, f"{sample_name}.bpNormalizedCoverage.Chr9.txt.gz")
    # save as txt
    df.to_csv(output_file, sep='\t', index=False, compression='gzip')
    
    print(f"Processed sample {sample_name} and saved the output in:  {output_file}")


if __name__ == "__main__":
    start_time = time.time()  # Record the start time
    
    # Arguments
    parser = argparse.ArgumentParser(description='Calculate normalized coverage at bp-resolution at give chromosome')
    parser.add_argument('-i', '--input', type=str, required=True, help='input file path')
    parser.add_argument('-c', '--chromosome', type=str, required=True, help='Name of target chromosome')
    args = parser.parse_args() 
    
    column_names = ['Chr', 'Position', 'DepthCoverage']
    sample_name = os.path.basename(args.input).replace("_Q30.depth.gz", "")
    
    # Check if input file exists
    check_file_exists(args.input)
    
    # Process the sample
    process_sample(args.input)
    
    end_time = time.time()
    duration = end_time - start_time
    
    print("Execution time:", duration, "seconds")
