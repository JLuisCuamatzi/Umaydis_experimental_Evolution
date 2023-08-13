# Genomic analysis in pools from line C

We sequenced pools along the experimental evolution, in order to have a lanscape of the main mutations detected in the colonies.

The genomes sequenced came from the next generations: 20, 30, 50, 70, 100, 140 and 200.

The `fastq` of the populations were cleaned and aligned to the reference genome with the same pipeline as the fastqs of the sequenced colonies.

We used the next python scripts: `01_Cleaning.py`, `02_Mapping.py`, and `02.1_NormalizedCoverage_and_Plotting.py`  to generate `sh` files to to process the fastqs
 - For this case, we used the next sample sheet: `USMA_EE_Pools_SampleSheet.csv`

<b>Cleaning</b>

```
# move to 04_Pools_Analysis

cd ~/Umaydis_experimental_Evolution/04_Pools_analysis/

python3 ../01_Cleaning/01_Cleaning.py -c ../   

# Execution of sh files:
for files in shFiles/*.sh; do sh $files; done

```

<b>Mapping And Coverage Plotting</b>

```

## Mapping

python3 ../02_MappingAndCoveragePlotting/02_Mapping.py -c ../USMA_EE_Pools_SampleSheet.csv



```