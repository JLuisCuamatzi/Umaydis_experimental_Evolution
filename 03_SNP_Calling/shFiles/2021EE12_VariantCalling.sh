#!/bin/bash

echo 'Start the variant calling for 2021EE12'

cd ~/Umaydis_experimental_Evolution/03_SNP_Calling/
# Defining reference genome
reference_genome="../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta"
# Defining inputs and outputs
input_cram_file="../02_MappingAndCoveragePlotting/bamFiles/2021EE12_BWA.mrkdup.addgp.cram" # File with aligment
output_mpileup_file="mpileupFiles/2021EE12_VC_bcftools.pileup"
output_vcf_file="vcfFiles/2021EE12_VC_bcftools.vcf.gz"
# create directories:
mkdir -p mpileupFiles/
mkdir -p vcfFiles/
# getting the mpileup file:
bcftools mpileup -f ${reference_genome} ${input_file} -q 20 -Q 30 > ${output_mpileup_file}
# getting vcf file:
bcftools call -m -v -Oz -o ${output_vcf_file} ${output_mpileup_file} --ploidy 1
# compressing pileup file:
gzip ${output_mpileup_file}

echo 'End the variant calling for 2021EE12'
#
