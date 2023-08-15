# Analysis of aneuploidies of <i>Ustilago maydis</i> isolates

We analyzed the ocurrence of chromosomal duplications in 50 <i>U. maydis</i> isolates reported in the next two studies:
 - [Population Genomics of the Maize Pathogen Ustilago maydis: Demographic History and Role of Virulence Clusters in Adaptation](https://academic.oup.com/gbe/article/13/5/evab073/6219951?login=true)
 - [Effectors with Different Gears: Divergence of Ustilago maydis Effector Genes Is Associated with Their Temporal Expression Pattern during Plant Infection](https://www.mdpi.com/2309-608X/7/1/16)
 
1.- The sequenced genomes were download from NCBI SRA

2.- These genomes were aligned to the reference genomes and the normalized coverage was obtained with the pipeline stated in `02_MappingAndCoverageTables`
 - The files with the normalized coverage are in: `~/Umaydis_experimental_Evolution/05_Duplications_in_Umaydis_isolates/normalizedCoverage/`

3.- The normalized coverage was plotted with the next r script

```
# The inputs for this script are the files in `~/Umaydis_experimental_Evolution/05_Duplications_in_Umaydis_isolates/normalizedCoverage/`

Rscript Duplications_Umaydis_isolates.R

```
