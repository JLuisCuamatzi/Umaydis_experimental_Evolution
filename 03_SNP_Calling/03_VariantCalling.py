"""
@ author: jcuamatzi

this script needs the next:
-c SampleSheet.csv


USAGE:

python3 03_VariantCalling.py -c ../USMA_EE_Colonies_SampleSheet.csv


"""

import csv
import os
import argparse


# Function to write the lines of sh file
def generate_sh_file(sampleID):
    # Generate the .sh file
    sh_path = "shFiles/"
    os.makedirs(sh_path, exist_ok=True)
    sh_name = f"{sampleID}_VariantCalling.sh"
    sh_file = os.path.join(sh_path, sh_name)
    #
    with open(sh_file, "w" ) as sh_file:
        sh_file.write(f'#!/bin/bash\n')
        sh_file.write(f'\n')
        sh_file.write(f"echo 'Start the variant calling for {sampleID}'\n")
        sh_file.write(f'\n')
        #
        
        sh_file.write(f'cd ~/Umaydis_experimental_Evolution/03_SNP_Calling/\n')
        #
        sh_file.write('# Defining reference genome\n')
        sh_file.write('reference_genome="../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta"\n')
        sh_file.write('# Defining inputs and outputs\n')
        sh_file.write(f'input_cram_file="../02_MappingAndCoveragePlotting/bamFiles/{sampleID}_BWA.mrkdup.addgp.cram" # File with aligment\n')
        sh_file.write(f'output_mpileup_file="mpileupFiles/{sampleID}_VC_bcftools.pileup"\n')
        sh_file.write(f'output_vcf_file="vcfFiles/{sampleID}_VC_bcftools.vcf.gz"\n')
        sh_file.write('# create directories:\n')
        sh_file.write('mkdir -p mpileupFiles/\n')
        sh_file.write('mkdir -p vcfFiles/\n')
        #
        sh_file.write('# getting the mpileup file:\n')
        sh_file.write('bcftools mpileup -f ${reference_genome} ${input_cram_file} -q 20 -Q 30 > ${output_mpileup_file}\n')
        #
        sh_file.write('# getting vcf file:\n')
        sh_file.write('bcftools call -m -v -Oz -o ${output_vcf_file} ${output_mpileup_file} --ploidy 1\n')
        #
        sh_file.write('# compressing pileup file:\n')
        sh_file.write('gzip ${output_mpileup_file}\n')
        sh_file.write('\n')
        #
        sh_file.write(f"echo 'End the variant calling for {sampleID}'\n")
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
    
    
    
