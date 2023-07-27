"""
@ author: jcuamatzi

this script needs the next:
-c SampleSheet.csv


USAGE:

python3 02.1_NormalizedCoverage_and_Plotting.py -c ../USMA_EE_Colonies_SampleSheet.csv


"""

import csv
import os
import argparse


# Function to write the lines of sh file
def generate_sh_file(sampleID):
    # Generate the .sh file
    sh_path = "shFiles/"
    os.makedirs(sh_path, exist_ok=True)
    sh_name = f"{sampleID}_NormCovAndPlotting.sh"
    sh_file = os.path.join(sh_path, sh_name)
    
    
    
    with open(sh_file, "w" ) as sh_file:
        sh_file.write(f'#!/bin/bash\n')
        sh_file.write(f'\n')
        sh_file.write(f"echo 'Start to Compute Normalized Coverage for {sampleID}'\n")
        #
        sh_file.write(f'# Calculate normalized coverage with the following python script: ComputeMedianCoverage.py\n')
        sh_file.write(f'FileDepthCoverage="coverageFiles/{sampleID}_Q30.depth.gz"\n')
        sh_file.write(f'FileNormalizedCoverage="normalizedCoverageTables/{sampleID}_NormalizedCoverage.txt"\n')
        sh_file.write(f'sample="{sampleID}"\n')
        #
        sh_file.write(f'# Compute Normalized Coverage with the next python script:ComputeMedianCoverage.py\n')
        sh_file.write('python3 ComputeMedianCoverage.py --input ${FileDepthCoverage} --output ${FileNormalizedCoverage} -w 1000\n')
        sh_file.write(f'#\n')
        #
        sh_file.write(f'# Create the coverage plots (raw and normalized) with the next R script: CoveragePlottR.R\n')
        sh_file.write('Rscript CoveragePlottR.R --normalizedCov_file ${FileNormalizedCoverage} --window_size 1000 --sample ${sample} --chr "USMA_521_v2_"\n')
        sh_file.write(f'\n')

        sh_file.write(f"#\n")


def main(csv_file):
    with open(csv_file, "r") as file:
        reader = csv.DictReader(file)
        for row in reader:
            sampleID = row["SampleID"]
            generate_sh_file(sampleID)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate .sh files")
    parser.add_argument("-c", "--csv_file", help="CSV file containing sample IDs", required=True) 
    args = parser.parse_args() 
    main(args.csv_file)
    
    
    
