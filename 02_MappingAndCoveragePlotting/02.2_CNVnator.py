"""
@ author: jcuamatzi

USAGE:

python3 02.2_CNVnator.py -c ../USMA_EE_Colonies_SampleSheet.csv

python3 02.2_CNVnator.py -c ../USMA_EE_Pools_SampleSheet.csv


"""

import csv
import os
import argparse


# Function to write the lines of sh file
def generate_sh_file(sampleID):
    # Generate the .sh file
    sh_path = "shFiles/"
    os.makedirs(sh_path, exist_ok=True)
    sh_name = f"{sampleID}_CNVnator.sh"
    sh_file = os.path.join(sh_path, sh_name)
    
    
    
    with open(sh_file, "w" ) as sh_file:
        sh_file.write(f'#!/bin/bash\n')
        sh_file.write(f'\n')
        sh_file.write(f"echo 'Start to Run CNVnator for {sampleID}'\n")
        #
        sh_file.write('# Create variables\n')
        sh_file.write("ref_genome='../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta'\n")
        sh_file.write(f"cramFile='bamFiles/{sampleID}_BWA.mrkdup.addgp.cram'\n")
        sh_file.write(f"rootFile='CNV_by_CNVnator/rootFiles/{sampleID}.root'\n")
        sh_file.write(f"outputFile='CNV_by_CNVnator/cnvFiles/{sampleID}.cnvnator'\n")
        #
        sh_file.write('# Create dirs\n')
        sh_file.write('mkdir -p CNV_by_CNVnator/rootFiles/\n')
        sh_file.write('mkdir -p CNV_by_CNVnator/cnvFiles/\n')
        #
        sh_file.write('cnvnator -root $rootFile -genome $ref_genome -tree $cramFile\n')
        sh_file.write('cnvnator -root $rootFile -his 150 -d ../../USMA_Genome/USMA_521_v2/\n')
        sh_file.write('cnvnator -root $rootFile -stat 150\n')
        sh_file.write('cnvnator -root $rootFile -partition 150\n')
        sh_file.write('cnvnator -root $rootFile -call 150 > $outputFile\n')
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
    
    
    
