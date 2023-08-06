## SNP Calling with `bcftools`

For the SNP calling, we used `bcftools`

The `cram` files with the reads from `fastq` aligned to <i>U. maydis</i> reference genome were the input for `bcftools mpileup`

The pipeline consist of three stages:

<b>03.1.- bcftools mpileup and call </b>

In this stage, we obtained a `pileup` file from the aligment (`cram` file) and subsequently, this `pileup` file was used to identifed variants

<b>03.2.- Remove SG200 Background in SNPs and Filter by QUAL </b>

In this stage, we filtered the SNPs by removing those that comes from the SG200 and then applying a quality filter

<b>03.3.- Filter SNPs by Depth Coverage in each sample </b>

Here, we work with the SNPs that were left as results of the previous step.
 - We check that, at the position where the SNP were identified, the alternative allele is supported by at least 90% of the reads aligned to that coordinate.
 - In addition, we removed those SNPs that in strain SG200, at least 40% of the aligned reads contained the identified alternative allele at that position.

<b>03.4.-  </b>




### 03.1.- bcftools mpileup and call

For this stage we used the next code:

```
# Running for sample: 2021EE01

echo 'Start mpileup and call for 2021EE01'

cd 03_SNP_Calling/

    # Defining reference genome
reference_genome="../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta"

    # Defining inputs and outputs
input_cram_file="../02_MappingAndCoveragePlotting/bamFiles/2021EE01_BWA.mrkdup.addgp.cram" # File with aligment
output_mpileup_file="mpileupFiles/2021EE01_VC_bcftools.pileup"
output_vcf_file="vcfFiles/2021EE01_VC_bcftools.vcf.gz"

    # create directories:
mkdir -p mpileupFiles/
mkdir -p vcfFiles/


# getting the mpileup file:
bcftools mpileup -f ${reference_genome} ${input_file} -q 20 -Q 30 > ${output_mpileup_file}

# getting vcf file:
bcftools call -m -v -Oz -o ${output_vcf_file} ${output_mpileup_file} --ploidy 1

# compressing 'pileup' file:
gzip ${output_mpileup_file}


```

To optimize the execution of this process, I wrote the `03_VariantCalling.py` python script that generates sh files with the aforementioned commands.

This script needs a sample sheet file in `csv` format where the ID of each sample is listed. 
 - The sample sheet is in the directory: `~/Umaydis_experimental_Evolution/USMA_EE_Colonies_SampleSheet.csv`.
 - `python3` is required to execute this script.

Execute `03_VariantCalling.py` script
```
# verify that you are in the directory: ~/Umaydis_experimental_Evolution/03_VariantCalling.py

python3 03_VariantCalling.py -c ../USMA_EE_Colonies_SampleSheet.csv

```

Once this script is executed, a directory named `shFiles` will be created. In this directory will be the `.sh` files to be executed.

You can run these `.sh` files with the next command:

```
# Example for sample 2021EE01

# with pwd, check that you are in: ~/Umaydis_experimental_Evolution/03_SNP_Calling, if not, move with cd to ~/Umaydis_experimental_Evolution/03_SNP_Calling
# pwd
# cd ~/Umaydis_experimental_Evolution/03_SNP_Calling

sh shFiles/2021EE01_VariantCalling.sh

```

The main output of this stage are individual `vcf` files for each analyzed sample. Also `mpileups` files were produced.

 - The `vcf` files are in `vcfFiles/` compressed as `.vcf.gz` files
 - The `pileup` files are in `mpileupFiles/` compressed as `.pileup.gz` files

 ### 03.2.- Remove SG200 Background in SNPs and Filter by QUAL
 
 The next step is filtering the SNPs as follows:
 - A) Identify and remove the SNPs already present in the departing strain <i>U. maydis</i> SG200 (using the data obtained after sequencing with Illumina and BGI)
 - B) Filter the remaining SNPs by QUAL (Q > 200)
 - C) Identify recurrent SNPs
 
To do these three steps we call the next R script: `SNP_Processing_and_Filtering.Step1.R`

Execute `SNP_Processing_and_Filtering.Step1.R` script <br>
 - This script requires the file: `Functions_Filtering_SNPs.R`, which is located in `~/Umaydis_experimental_Evolution/03_SNP_Calling/`

```
# Verify that you are in the directory: ~/Umaydis_experimental_Evolution/03_SNP_Calling/

Rscript SNP_Processing_and_Filtering.Step1.R

```
 
After the execution of this script, the ouptus (tables & plots) were writted in:
 - `~/Umaydis_experimental_Evolution/03_SNP_Calling/SNPCalling`<br>
     `|-Plots/`<br>
     `|-Tables/`<br>

### 03.3.- Filter SNPs by Depth Coverage in each sample

Once we have already filtered the SNPs by removing the SG200 background and by removing those with a QUAL minor than 200, the next step is to check if at the given SNP position, such SNP is supported by at least 90% of the reads aligned to that position.

For this purpose, we extracted the information of reads aligned at the SNP position using `bcftools` with the next code:

```
# Example for sample: 2021EE01

# Indicate the bed file
bedFile="SNPCalling/Tables/SNPsToCheck.bed"        

# Indicate the bam file (in this case compressed as cram)
bamFile="../02_MappingAndCoveragePlotting/bamFiles/2021EE01_BWA.mrkdup.addgp.cram"

# Indicate the reference genome
reference_genome="../USMA_Genome/USMA_521_v2/USMA_521_v2.fasta"

# Create the directory for the output and indicate the output file
mkdir -p SNPCalling/Tables/Check_SNPs_AF
outputFile=SNPCalling/Tables/Check_SNPs_AF/2021EE01_SNPs_AlleleDepth.txt

# Extract the aligned nucleotides at the position indicated in the bed file
bcftools mpileup -q 20 -Q 30 --regions-file ${bed_file} ${bam_file} --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --fasta-ref ${reference_genome} | bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > ${outputFile}

echo 'End the checking SNPs AF in 2021EE01'

```

To optimize the execution of this process, I wrote the `03.1_Check_SNPs_AF.py` python script that generates sh files with the aforementioned commands.

Execute the python script
```
python3 03.1_Check_SNPs_AF.py -c ../USMA_EE_Colonies_SampleSheet.csv

```

Once this script is executed, the `.sh` files will be writtend in `shFiles/`.

You can run these `.sh` files with the next command:

```
# Example for sample 2021EE01

# with pwd, check that you are in: ~/Umaydis_experimental_Evolution/03_SNP_Calling, if not, move with cd to ~/Umaydis_experimental_Evolution/03_SNP_Calling

sh shFiles/2021EE01_Check_SNPs_AF.sh

```

The output of this step are `txt` files with the counting of each nucleotide aligned at this position.
 - These files are in: `SNPCalling/Tables/Check_SNPs_AF`
 
Now, we filtered the SNPs with the next criteria:
 - Keep those SNPs where the alternative allele is supported by at least 90% of the reads aligned to that coordinate.
 - Remove those SNPs that in strain SG200, at least 40% of the aligned reads contained the identified alternative allele at that position.
 
For this task, we used the next custom R script: `SNP_Processing_and_Filtering.Step2.R`
 This script requires the next files produced in `SNP_Processing_and_Filtering.Step1.R` as input:
   - `Shared.SNP.NoSG200.Q200.csv`: matrix of shared SNPs in colonies among lines
   - `Filtered.Variants.Q200.csv`: data set wtih SNPs 
   - `SNP_Number_Per_Sample.Q200.csv`: data set that summarizes the number of SNPs remaining after each filtering step
   
```
# Execution of: SNP_Processing_and_Filtering.Step2.R
Rscript SNP_Processing_and_Filtering.Step2.R

```

The outputs of this script are:
 The next tables located in `~/Umaydis_experimental_Evolution/03_SNP_Calling/SNPCalling/Tables/`
  - `Shared.SNP.NoSG200.Q200.AF90.csv`
  - `SNP_Number_Per_Sample.Q200.AF90.csv`
  - `Matriz.Q200.AF90.Raw.csv`
  - `Matriz.Q200.AF90.Percentage.csv`
  <br>
 And the next plots located in: `~/Umaydis_experimental_Evolution/03_SNP_Calling/SNPCalling/Plots/`
  - `Matrix.Q200.AF90.RAW.svg`
  - `Matrix.Q200.AF90.PER.svg`
 
 <b> The final list of SNPs is in `Shared.SNP.NoSG200.Q200.AF90.csv` </b>
 
 
 










