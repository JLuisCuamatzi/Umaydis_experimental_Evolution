#!/bin/bash

echo 'Start cleaning for 2021EE22'

# Inputs
input_fastqR1="fastqFiles/raw/2021EE22_R1.fq.gz"
input_fastqR2="fastqFiles/raw/2021EE22_R2.fq.gz"

# Outputs
output_fastqR1="fastqFiles/clean/2021EE22_R1_clean.fastq.gz"
output_fastqR2="fastqFiles/clean/2021EE22_R2_clean.fastq.gz"

# Unpaired fastq files
unpaired1="fastqFiles/unpaired/2021EE22_R1_unpaired.fastq.gz"
unpaired1="fastqFiles/unpaired/2021EE22_R2_unpaired.fastq.gz"

# json file
json_file="fastpReports/2021EE22_fastp.json" 

# html file
html_file="fastpReports/2021EE22_fastp.html" 

fastp -i ${input_fastqR1} -o ${output_fastqR1} -I ${input_fastqR2} -O ${output_fastqR2} --unpaired1 ${unpaired1} --unpaired2 ${unpaired2} -w 2 -y -x -z 9 -j ${json_file} -h ${html_file} 

echo 'End cleaning for 2021EE22'
