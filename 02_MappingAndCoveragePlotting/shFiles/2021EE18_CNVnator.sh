#!/bin/bash

echo 'Start to Run CNVnator for 2021EE18'
# Create variables
ref_genome='../../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta'
cramFile='../bamFiles/2021EE10_BWA.mrkdup.addgp.cram'
rootFile='rootFiles/2021EE10.root'
outputFile='cnvFiles/2021EE10.cnvnator'
# Create dirs
mkdir -p rootFiles/
mkdir -p cnvFiles/
cnvnator -root $rootFile -genome $ref_genome -tree $cramFile
cnvnator -root $rootFile -his 150 -d ../../USMA_Genome/USMA_521_v2/
cnvnator -root $rootFile -stat 150
cnvnator -root $rootFile -partition 150
cnvnator -root $rootFile -call 150 > $outputFile
#
