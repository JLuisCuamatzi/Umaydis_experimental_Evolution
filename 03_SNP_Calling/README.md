## SNP Calling with `bcftools`

For the SNP calling, we used `bcftools`

The `cram` files with the reads from `fastq` aligned to <i>U. maydis</i> reference genome were the input for `bcftools mpileup`

The pipeline consist of three stages:

<b>03.1.- bcftools mpileup and call </b>

In this stage we obtained a `pileup` file from the aligment (`cram` file) and subsequently, this `pileup` file was used to identifed variants

<b>03.2.- Remove SG200 Background in SNPs and Filter by QUAL </b>



<b>03.3.- Filter SNPs by Depth Coverage in each sample </b>






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
python3 03_VariantCalling.py -c ../USMA_EE_Colonies_SampleSheet.csv

```

Once this script is executed, a directory named `shFiles` will be created. In this directory will be the `.sh` files to be executed.

You can run these `.sh` files with the next command:

```
# Example for sample 2021EE01

sh shFiles/2021EE01_VariantCalling.sh

```

The main output of this stage are individual VCF for each analyzed sample. Also mpileups files were produced.

 - The `vcf` files are in `vcfFiles/` compressed as `.vcf.gz` files
 - The `pileup` files are in `mpileupFiles/` compressed as `.pileup.gz` files

 ### 03.2.- Remove SG200 Background in SNPs and Filter by QUAL
 
 The next step is filtering the SNPs as follows:
 - A) Identify and remove the SNPs already present in the departing strain <i>U. maydis</i> SG200 (using the data obtained after sequencing with Illumina and BGI)
 - B) Filter the remaining SNPs by QUAL (Q > 200)
 - C) Identify recurrent SNPs
 
To do these three steps we call the next R script: `SNP_Processing_and_Filtering.Step1.R`

Execute `SNP_Processing_and_Filtering.Step1.R` script
 - This script requires the file: `Functions_Filtering_SNPs.R`

```

Rscript --vanilla SNP_Processing_and_Filtering.Step1.R --workingDir /mnt/Guanina/lmorales/Public/Ustilago_experimental_evolution_2021/

```
 











