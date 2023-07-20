"""
@ author: jcuamatzi

this script needs the next:
-c SampleSheet.csv


USAGE:

python3 02_Cleaning.py -c ../USMA_EE_Colonies_SampleSheet.csv


"""

import csv
import os
import argparse


# Function to write the lines of sh file
def generate_sh_file(sampleID):
    # Generate the .sh file
    sh_path = "shFiles/"
    os.makedirs(sh_path, exist_ok=True)
    sh_name = f"{sampleID}_Mapping.sh"
    sh_file = os.path.join(sh_path, sh_name)
    # Bam directories
    bam_path = "bamFiles"
    os.makedirs(bam_path, exist_ok=True)
    bam_stats_path = "bamFiles/stats/"
    os.makedirs(bam_stats_path, exist_ok=True)
    coverage_path = "coverageFiles"
    os.makedirs(coverage_path, exist_ok=True)
    
    
    with open(sh_file, "w" ) as sh_file:
        sh_file.write(f'#!/bin/bash\n')
        sh_file.write(f'\n')
        sh_file.write(f"echo 'Start mapping for {sampleID}'\n")
        #
        sh_file.write(f'\n')
        sh_file.write(f'# Defining reference genome\n')
        sh_file.write(f'reference_genome="../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta"\n')
        sh_file.write(f'\n')
        #
        sh_file.write(f'# Defining inputs and outputs\n')
        sh_file.write(f'fq_file_1="../01_Cleaning/fastqFiles/clean/{sampleID}_R1_clean.fastq.gz"\n') # fq R1
        sh_file.write(f'fq_file_2="../01_Cleaning/fastqFiles/clean/{sampleID}_R2_clean.fastq.gz"\n') # fq R2
        #
        sh_file.write(f'\n')
        sh_file.write(f'bam_file="bamFiles/{sampleID}_BWA.bam"\n')  # bam file
        #
        sh_file.write(f'\n')
        sh_file.write(f'bam_mrkdup_file="bamFiles/{sampleID}_BWA.mrkdup.bam"\n')  # bam file mark duplicates
        sh_file.write(f'\n')
        sh_file.write(f'dupMtrx_file="bamFiles/stats/{sampleID}_BWA_DuplicateMatrix"\n') # metrics file
        #
        sh_file.write(f'\n')
        sh_file.write(f'bam_mrkdup_addgp_file="bamFiles/{sampleID}_BWA.mrkdup.addgp.bam"\n')  # bam file
        # in cram
        sh_file.write(f'\n')
        sh_file.write(f'cram_mrkdup_addgp_file="bamFiles/{sampleID}_BWA.mrkdup.addgp.cram"\n')  # cram file 
        # depth file
        sh_file.write(f'\n')
        sh_file.write(f'depth_file="coverageFiles/{sampleID}_Q30.depth"\n')
        sh_file.write(f'\n')       
        sh_file.write(f"echo 'Running mapping with bwa mem for {sampleID}'\n")
        sh_file.write(f'\n')
        sh_file.write('''bwa mem -M -t10 ${reference_genome} ${fq_file_1} ${fq_file_2} | samtools view -hbS - | samtools sort -o ${bam_file} - \n''')
        sh_file.write(f'\n')
        sh_file.write(f'# Mark duplicates with picard\n')
        sh_file.write('''picard MarkDuplicates INPUT=${bam_file} OUTPUT=${bam_mrkdup_file} METRICS_FILE=${dupMtrx_file} VALIDATION_STRINGENCY=LENIENT\n''')
        sh_file.write(f'\n')
        #
        sh_file.write(f'# Add groups\n')
        sh_file.write('picard AddOrReplaceReadGroups I=${bam_mrkdup_file} O=${bam_mrkdup_addgp_file} LB=${sampleID} PL=illumina PU=${sampleID} SM=${sampleID} VALIDATION_STRINGENCY=LENIENT\n')
        sh_file.write(f'\n')
        sh_file.write('samtools index ${bam_mrkdup_addgp_file}\n')
        sh_file.write(f'\n')
        sh_file.write(f'echo "End mapping for {sampleID}"\n')
        # pass from bam to cram
        sh_file.write(f'\n')
        sh_file.write(f'# Compress bam to cram\n')
        sh_file.write('samtools view -T ${reference_genome} -C -o ${cram_mrkdup_addgp_file} ${bam_mrkdup_addgp_file}\n')
        sh_file.write('samtools index ${cram_mrkdup_addgp_file}\n')
        # safely remove bam files
        sh_file.write(f'\n')
        sh_file.write(f'# Remove bam files\n')
        sh_file.write('if [[ -s ${cram_mrkdup_addgp_file} ]]; then rm -rf ${bam_file}; fi\n')
        sh_file.write('if [[ -s ${cram_mrkdup_addgp_file} ]]; then rm -rf ${bam_mrkdup_file}; fi\n')
        sh_file.write('if [[ -s ${cram_mrkdup_addgp_file} ]]; then rm -rf ${bam_mrkdup_addgp_file}; fi\n')
        # Extract Coverage
        sh_file.write(f'# Extract coverage from alignment using samtools\n')
        sh_file.write('samtools depth -a -Q 30 ${cram_mrkdup_addgp_file} > ${depth_file}\n') # reads with Q > 30
        sh_file.write('gzip ${depth_file}\n') # compress `.depth` file
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
    
    
    
