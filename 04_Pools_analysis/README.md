# Genomic analysis in pools from line C

We sequenced pools along the experimental evolution, in order to have a lanscape of the main mutations detected in the colonies.

The genomes sequenced came from the next generations: 20, 30, 50, 70, 100, 140 and 200.

On these genomes, we looked for:
 1.- SNP in UMAG_05545 (Chr18:330,272)
 2.- Chromosome 9 duplication


First, we created a `bed` file that containts the next coordinates: Chr18:330,272. For this we can do a subset from the bed file at: `../03_SNP_Calling/SNPCalling/Tables/SNPsToCheck.bed`

```
# in bash terminal:

head -n 1 ../03_SNP_Calling/SNPCalling/Tables/SNPsToCheck.bed > SNP_in_UMAG_05545.bed

```

The file `SNP_in_UMAG_05545.bed` is the input `bcftools mpileup` as follows:


```
reference_genome="../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta"
bed_file="SNP_in_UMAG_05545.bed"

# Loop with bcftools mpileup for each sample

for sample_id in 2021EE30 2021EE31 2021EE32 2021EE33 2021EE34 2021EE35 2021EE36; do
    bam_file="../02_MappingAndCoveragePlotting/bamFiles/${sample_id}_BWA.mrkdup.addgp.cram"
    outputFile="${sample_id}_SNP_UMAG_05545.txt"

    bcftools mpileup -q 20 -Q 30 --regions-file ${bed_file} ${bam_file} --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --fasta-ref ${reference_genome} | bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > ${outputFile}
done

```

Then the Duplication of Chr 9 and the SNP in UMAG_05545 is analyzed with the next script: `CheckingVariantsInPools.R`
This script requires as input:
 - Files with the normalized coverage in the pools. These files are in: `../02_MappingAndCoveragePlotting/normalizedCoverageTables/`
 - Files with the count of alleles at the position chr18:330272. These files are in: `04_Pools_analysis`

```

Rscript CheckingVariantsInPools.R

```

The output of this file is the next file: `Figure_Variants.Pools.png`
