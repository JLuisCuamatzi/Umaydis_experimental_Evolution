"""
@ author: jcuamatzi

this script needs the next:
-c SampleSheet.csv


USAGE:

python3 03.1_Check_SNPs_AF.py -c ../USMA_EE_Colonies_SampleSheet.csv


"""

import csv
import os
import argparse

# Function to write the lines of sh file
def generate_sh_file(sampleID):
    # Generate the .sh file
    sh_path = "shFiles/"
    os.makedirs(sh_path, exist_ok=True)
    sh_name = f"{sampleID}_Check_SNPs_AF.sh"
    sh_file = os.path.join(sh_path, sh_name)
    #
    with open(sh_file, "w" ) as sh_file:
        sh_file.write(f'#!/bin/bash\n')
        sh_file.write(f'\n')
        sh_file.write(f"echo 'Start the checking the AF in {sampleID}'\n")
        sh_file.write(f'\n')
        # bed file        
        sh_file.write('bedFile="SNPCalling/Tables/SNPsToCheck.bed"\n')
        sh_file.write(f'\n')
        # bam file
        sh_file.write(f'bamFile="../02_MappingAndCoveragePlotting/bamFiles/{sampleID}_BWA.mrkdup.addgp.cram"\n')
        sh_file.write(f'\n')
        # reference genome
        sh_file.write('reference_genome="../../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta"\n')
        sh_file.write(f'\n')
        # output file
        sh_file.write('mkdir -p SNPCalling/Tables/Check_SNPs_AF\n')
        sh_file.write(f'outputFile=SNPCalling/Tables/Check_SNPs_AF/{sampleID}_SNPs_AlleleDepth.txt\n')
        sh_file.write(f'\n')
        #   
        sh_file.write('''bcftools mpileup -q 20 -Q 30 --regions-file ${bed_file} ${bam_file} --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --fasta-ref ${reference_genome} | bcftools query --format '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n' > ${outputFile} \n''')
        sh_file.write('\n')
        #
        sh_file.write(f"echo 'End the checking SNPs AF in {sampleID}'\n")
        sh_file.write("#\n")


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
    
    
