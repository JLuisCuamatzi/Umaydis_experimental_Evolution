#!/bin/bash

echo 'Start mapping for 2021EE36'

# Defining reference genome
reference_genome="../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta"

# Defining inputs and outputs
fq_file_1="../01_Cleaning/fastqFiles/clean/2021EE36_R1_clean.fastq.gz"
fq_file_2="../01_Cleaning/fastqFiles/clean/2021EE36_R2_clean.fastq.gz"

bam_file="bamFiles/2021EE36_BWA.bam"

bam_mrkdup_file="bamFiles/2021EE36_BWA.mrkdup.bam"

dupMtrx_file="bamFiles/stats/2021EE36_BWA_DuplicateMatrix"

bam_mrkdup_addgp_file="bamFiles/2021EE36_BWA.mrkdup.addgp.bam"

cram_mrkdup_addgp_file="bamFiles/2021EE36_BWA.mrkdup.addgp.cram"

depth_file="coverageFiles/2021EE36_Q30.depth"

echo 'Running mapping with bwa mem for 2021EE36'

bwa mem -M -t10 ${reference_genome} ${fq_file_1} ${fq_file_2} | samtools view -hbS - | samtools sort -o ${bam_file} - 

# Mark duplicates with picard
picard MarkDuplicates INPUT=${bam_file} OUTPUT=${bam_mrkdup_file} METRICS_FILE=${dupMtrx_file} VALIDATION_STRINGENCY=LENIENT

# Add groups
picard AddOrReplaceReadGroups I=${bam_mrkdup_file} O=${bam_mrkdup_addgp_file} LB=${sampleID} PL=illumina PU=${sampleID} SM=${sampleID} VALIDATION_STRINGENCY=LENIENT

samtools index ${bam_mrkdup_addgp_file}

echo "End mapping for 2021EE36"

# Compress bam to cram
samtools view -T ${reference_genome} -C -o ${cram_mrkdup_addgp_file} ${bam_mrkdup_addgp_file}
samtools index ${cram_mrkdup_addgp_file}

# Remove bam files
if [[ -s ${cram_mrkdup_addgp_file} ]]; then rm -rf ${bam_file}; fi
if [[ -s ${cram_mrkdup_addgp_file} ]]; then rm -rf ${bam_mrkdup_file}; fi
if [[ -s ${cram_mrkdup_addgp_file} ]]; then rm -rf ${bam_mrkdup_addgp_file}; fi
# Extract coverage from alignment using samtools
samtools depth -a -Q 30 ${cram_mrkdup_addgp_file} > ${depth_file}
gzip ${depth_file}
#
