## 02.- Mapping and Coverage Plotting

In the next stage of our pipeline, we did the following:

02.1.- We aligned the short-reads to the reference genome of <i>U. maydis</i>, and subsequently, we extracted the depth coverage of each base.

02.2.- We computed and plotted the normalized coverage of each genome in non-overlapping windows of 1 kb.


### 02.1.- Mapping

In this stept, the reads in `fastq` files were aligned to the reference genome of <i>Ustilago maydis</i>

For the alignment we used [BWA_mem]([https://academic.oup.com/bioinformatics/article/26/5/589/211735) using the next code:

```
## Running for sample: 2021EE01

echo 'Start mapping for 2021EE01'


cd 02_MappingAndCoveragePlotting/
    # Defining reference genome
reference_genome="../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta"
    # Defining inputs and outputs
fq_file_1="../01_Cleaning/fastqFiles/clean/2021EE01_R1_clean.fastq.gz"
fq_file_2="../01_Cleaning/fastqFiles/clean/2021EE01_R2_clean.fastq.gz"
bam_file="bamFiles/2021EE01_BWA.bam"
bam_mrkdup_file="bamFiles/2021EE01_BWA.mrkdup.bam"
dupMtrx_file="bamFiles/stats/2021EE01_BWA_DuplicateMatrix"
bam_mrkdup_addgp_file="bamFiles/2021EE01_BWA.mrkdup.addgp.bam"
cram_mrkdup_addgp_file="bamFiles/2021EE01_BWA.mrkdup.addgp.cram"
depth_file="coverageFiles/2021EE01_Q30.depth"


echo 'Running mapping with bwa mem for 2021EE01'
bwa mem -M -t10 ${reference_genome} ${fq_file_1} ${fq_file_2} | samtools view -hbS - | samtools sort -o ${bam_file} - 

# Mark duplicates with picard
picard MarkDuplicates INPUT=${bam_file} OUTPUT=${bam_mrkdup_file} METRICS_FILE=${dupMtrx_file} VALIDATION_STRINGENCY=LENIENT

# Add groups
picard AddOrReplaceReadGroups I=${bam_mrkdup_file} O=${bam_mrkdup_addgp_file} LB=${sampleID} PL=illumina PU=${sampleID} SM=${sampleID} VALIDATION_STRINGENCY=LENIENT

samtools index ${bam_mrkdup_addgp_file}

echo "End mapping for 2021EE01"

# Compress bam to cram
samtools view -T ${reference_genome} -C -o ${cram_mrkdup_addgp_file} ${bam_mrkdup_addgp_file}
samtools index ${cram_mrkdup_addgp_file}

# Remove bam files if cram already exists
if [[ -s ${cram_mrkdup_addgp_file} ]]; then rm -rf ${bam_file}; fi
if [[ -s ${cram_mrkdup_addgp_file} ]]; then rm -rf ${bam_mrkdup_file}; fi
if [[ -s ${cram_mrkdup_addgp_file} ]]; then rm -rf ${bam_mrkdup_addgp_file}; fi

# Extract coverage from alignment using samtools
samtools depth -a -Q 30 ${cram_mrkdup_addgp_file} > ${depth_file}
gzip ${depth_file}


```


To optimize the execution of this process, I wrote the `02_Mapping.py` python script that generates sh files with the aforementioned commands.

This script needs a sample sheet file in `csv` format where the ID of each sample is listed. 
 - The sample sheet is in the main directory: `~/Umaydis_experimental_Evolution/USMA_EE_Colonies_SampleSheet.csv`.
 - `python3` is required to execute this script.

Execute `02_Mapping.py` script
```
python3 01_Mapping.py -c ../USMA_EE_Colonies_SampleSheet.csv

```

Once this script is executed, a directory named `shFiles` will be created. In this directory will be the `.sh` files to be executed.

You can run these `.sh` files with the next command:

```
# Example for sample 2021EE01

sh shFiles/2021EE01_Mapping.sh

```

The output files of this process are:
  - `cram` files located in `bamFiles/`
  - `depth.gz` files located in `coverageFiles/` with the depth of coverage of each base.

### 02.2.- Computing and Plotting Normalized Coverage

The next step is calculate the normalized coverage in non-overlapping windows of 1 kb.
   - Normalized coverage: median coverage depth of each windows / global median coverage depth
   
 For this task, we use a <b>python script</b> to compute the normalized coverage and a <b>R script</b> to plot the raw and normalized coverage
 
```
echo 'Start to Compute Normalized Coverage for 2021EE01'

# variables:
FileDepthCoverage="coverageFiles/2021EE01_Q30.depth.gz"
FileNormalizedCoverage="normalizedCoverageTables/2021EE01_NormalizedCoverage.txt"
sample="2021EE01"


# Compute Normalized Coverage with the next python script: ComputeMedianCoverage.py
python3 ComputeMedianCoverage.py --input ${FileDepthCoverage} --output ${FileNormalizedCoverage} -w 1000


# Create the coverage plots (raw and normalized) with the next R script: CoveragePlottR.R
Rscript CoveragePlottR.R --normalizedCov_file ${FileNormalizedCoverage} --window_size 1000 --sample ${sample} --chr "USMA_521_v2_"

```

To optimize the execution of this process, I wrote the `02.1_NormalizedCoverage_and_Plotting.py` python script that generates sh files with the aforementioned commands.

Execute `02.1_NormalizedCoverage_and_Plotting.py` script

```
python3 02.1_NormalizedCoverage_and_Plotting.py -c ../USMA_EE_Colonies_SampleSheet.csv

```

Once this script is executed, the `.sh` files will be saved in `shFiles`.

You can run these `.sh` files with the next command:

```
# Example for sample 2021EE01

sh shFiles/2021EE01_NormCovAndPlotting.sh

```

### Ratio in coverage 

Finally, the ratio between coverage in each sample against reference sample (<i>U. maydis</i> SG200) was computed in R in the next script:

`Coverage_Ratio_Log2.Rmd`









