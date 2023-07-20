"""
@ author: jcuamatzi

this script needs the next:
-c SampleSheet.csv


USAGE:

python3 01_Cleaning.py -c ../USMA_EE_Colonies_SampleSheet.csv


"""

import csv
import os
import argparse


# Function to write the lines of sh file
def generate_sh_file(sampleID):
    # Generate the .sh file
    sh_path = "shFiles/"
    os.makedirs(sh_path, exist_ok=True)
    sh_name = f"{sampleID}_Cleaning.sh"
    sh_file = os.path.join(sh_path, sh_name)
    # Fastq Files
    raw_path = "fastqFiles/raw/"
    os.makedirs(raw_path, exist_ok=True)
    clean_path = "fastqFiles/clean/"
    os.makedirs(clean_path, exist_ok=True)
    unpaired_path = "fastqFiles/unpaired"
    os.makedirs(unpaired_path, exist_ok=True)
    report_path = "fastpReports/"
    os.makedirs(report_path, exist_ok=True)
    
    with open(sh_file, "w" ) as sh_file:
        sh_file.write(f'#!/bin/bash\n')
        sh_file.write(f'\n')
        sh_file.write(f"echo 'Start cleaning for {sampleID}'\n")
        #
        sh_file.write(f'\n')
        sh_file.write(f'# Inputs\n')
        sh_file.write(f'input_fastqR1="fastqFiles/raw/{sampleID}_R1.fq.gz"\n')
        sh_file.write(f'input_fastqR2="fastqFiles/raw/{sampleID}_R2.fq.gz"\n')
        #
        sh_file.write(f'\n')
        sh_file.write(f'# Outputs\n')
        sh_file.write(f'output_fastqR1="fastqFiles/clean/{sampleID}_R1_clean.fastq.gz"\n')
        sh_file.write(f'output_fastqR2="fastqFiles/clean/{sampleID}_R2_clean.fastq.gz"\n')
        #
        sh_file.write(f'\n')
        sh_file.write(f'# Unpaired fastq files\n')
        sh_file.write(f'unpaired1="fastqFiles/unpaired/{sampleID}_R1_unpaired.fastq.gz"\n')
        sh_file.write(f'unpaired1="fastqFiles/unpaired/{sampleID}_R2_unpaired.fastq.gz"\n')
        #
        sh_file.write(f'\n')
        sh_file.write(f'# json file\n')
        sh_file.write(f'json_file="fastpReports/{sampleID}_fastp.json" \n')
        #
        sh_file.write(f'\n')
        sh_file.write(f'# html file\n')
        sh_file.write(f'html_file="fastpReports/{sampleID}_fastp.html" \n')
        #
        sh_file.write(f'\n')
        sh_file.write('''fastp -i ${input_fastqR1} -o ${output_fastqR1} -I ${input_fastqR2} -O ${output_fastqR2} --unpaired1 ${unpaired1} --unpaired2 ${unpaired2} -w 2 -y -x -z 9 -j ${json_file} -h ${html_file} \n''')
        #
        sh_file.write(f'\n')
        sh_file.write(f"echo 'End cleaning for {sampleID}'\n")


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
    
    
    
