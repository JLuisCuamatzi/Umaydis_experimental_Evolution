#!/bin/bash

echo 'Start the checking the AF in 2021EE22'

bedFile="SNPCalling/Tables/SNPsToCheck.bed"

bamFile="../02_MappingAndCoveragePlotting/bamFiles/2021EE22_BWA.mrkdup.addgp.cram"

reference_genome="../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta"

mkdir -p SNPCalling/Tables/Check_SNPs_AF
outputFile=SNPCalling/Tables/Check_SNPs_AF/2021EE22_SNPs_AlleleDepth.txt

bcftools mpileup -q 20 -Q 30 --regions-file ${bed_file} ${bam_file} --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --fasta-ref ${reference_genome} | bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > ${outputFile} 

echo 'End the checking SNPs AF in 2021EE22'
#
