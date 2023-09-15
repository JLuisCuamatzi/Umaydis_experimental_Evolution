## 02.- Mapping and Coverage Plotting

In the next stage of our pipeline, we did the following:

<b>02.1.- Mapping</b>
 - We aligned the short-reads to the reference genome of <i>U. maydis</i>, and subsequently, from the alignment we extracted the depth coverage of each base.

<b>02.2.- Computing and Plotting Normalized Coverage </b>
 - We computed and plotted the normalized coverage of each genome in non-overlapping windows of 1 kb.
 
<b>02.3.- Computing the Ratio of Normalized Coverage </b>
 - We computed the ratio between normalized coverage in each sample against the initial strain (<i>U. maydis</i> SG200)
 
<b>02.4.- Mapping genomes from sequenced pools and analysis of coverage</b>

<b>02.5.- Analysis of coverage with a base-pair resolution at the breakpoint of chr 9 duplication</b>

<b>02.6.- Analysis of Copy Number Variations with CNVnator</b>

### 02.1.- Mapping

In this stept, the reads in `fastq` files were aligned to the reference genome of <i>Ustilago maydis</i>

For the alignment we used <b>[BWA_mem]([https://academic.oup.com/bioinformatics/article/26/5/589/211735)</b> using the next code:

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
python3 02_Mapping.py -c ../USMA_EE_Colonies_SampleSheet.csv

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

### 02.3.- Computing the Ratio of Normalized Coverage

Finally, the ratio between normalized coverage in each sample against reference sample (<i>U. maydis</i> SG200) was computed in R in the next script:

`Coverage_Ratio_Log2.R`

The file `02_MappingAndCoveragePlotting/Figures_RatioLog2/Plot_USMA_521_v2_9.Log2Ratio.png` is the <b>Figure 3</b> in the main manuscript.



### 02.4.- Mapping genomes from sequenced pools and analysis of coverage

The `fastq` of the pools were aligned with the same pipeline

We used the next python scripts: `02_Mapping.py` to generate `sh` files to to process these fastqs
 - For this case, we used the next sample sheet: `USMA_EE_Pools_SampleSheet.csv`

<b>Cleaning and plotting normalized coverage</b>

```
# To map the reads

python3 02_Mapping.py -c ../USMA_EE_Pools_SampleSheet.csv

# To calculate and plot the normalized coverage

python3 02.1_NormalizedCoverage_and_Plotting.py -c ../USMA_EE_Pools_SampleSheet.csv

# the sh files from both python scripts are in `shFiles/`, and can be executed with:

sh shFiles/2021EE30_Mapping.sh              # example for 2021EE30
sh shFiles/2021EE30_NormCovAndPlotting.sh   # example for 2021EE30

```

## 02.5.- Analysis of coverage with a base-pair resolution at the breakpoint of chr 9 duplication

An analysis of normalized coverage at base-pair resolution (bpRes) was done to identify the breakpoint of aneuploidy in chromosome 9

The normalized coverage at bpRes was computed from depth coverage files using the next python script: `ComputeMedianCoverage_bpResolution.py`
 - The input for this script is a `.depth.gz` file. The input should be indicated with `-i`
 - For computational efficiency, the analysis is just performed on a target chromosome. This should be indicated with `-c`
 

```
# Check in the terminal that the current directory is: ~/Umaydis_experimental_Evolution/02_MappingAndCoveragePlotting/

# execute the python for all depth coverage files (located in: coverageFiles/)

for files in coverageFiles/*.depth.gz; do python3 ComputeMedianCoverage_bpResolution.py -i $files -c USMA_521_v2; done


```

The output files are in: `normalizedCoverageTables/bpResolution/`

The normalized coverage at bp Resolution was plotted with the next Rscript: `NormalizedCoverage_bpRes.R`

Additionally, we plotted the genes flankig the breakpoint, and the GC percent in that region.


## 02.6.- Analysis of Copy Number Variations with CNVnator

We conducted an analysis with <b>[CNVnator](https://genome.cshlp.org/content/21/6/974.short)</b> in order to get a more precise size of the chromosomal amplifications detected by analyzing the log<sub>2</sub> ratio

First, is required to generate individual files for each chr in the reference genome

For this, we used the next Python script: `Split.Reference.Genome.py`

```
python3 Split.Reference.Genome.py -f ../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta

```

## Running CNVnator

The version used was `cnvator v.0.3.3`

We used the next code to process one sample:

```
# Create variables
# example for sample 2021EE01

ref_genome='../../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta'
cramFile='../bamFiles/2021EE01_BWA.mrkdup.addgp.cram'
rootFile='rootFiles/2021EE01.root'
outputFile='cnvFiles/2021EE01.cnvnator'

# create dirs
mkdir -p rootFiles/
mkdir -p cnvFiles/

# Execute CNVnator
cnvnator -root $rootFile -genome $ref_genome -tree $cramFile
cnvnator -root $rootFile -his 150 -d ../../USMA_Genome/USMA_521_v2/
cnvnator -root $rootFile -stat 150
cnvnator -root $rootFile -partition 150
cnvnator -root $rootFile -call 150 > $outputFile

```

To optimize the execution of this process, I wrote the `02.2_CNVnator.py` python script that generates sh files with the aforementioned commands.

```
python3 02.2_CNVnator.py -c ../USMA_EE_Colonies_SampleSheet.csv

python3 02.2_CNVnator.py -c ../USMA_EE_Pools_SampleSheet.csv

```

The sh files resulting from that command are located in `shFiles/`

Theses files can be executed with:

```
# Example for sample 2021EE01

sh shFiles/2021EE01_CNVnator.sh

```
