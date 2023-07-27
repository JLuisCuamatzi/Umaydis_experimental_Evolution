#!/bin/bash

echo 'Start to Compute Normalized Coverage for 2021EE12'
# Calculate normalized coverage with the following python script: ComputeMedianCoverage.py
FileDepthCoverage="coverageFiles/2021EE12_Q30.depth.gz"
FileNormalizedCoverage="normalizedCoverageTables/2021EE12_NormalizedCoverage.txt"
sample="2021EE12"
# Compute Normalized Coverage with the next python script:ComputeMedianCoverage.py
python3 ComputeMedianCoverage.py --input ${FileDepthCoverage} --output ${FileNormalizedCoverage} -w 1000
#
# Create the coverage plots (raw and normalized) with the next R script: CoveragePlottR.R
Rscript CoveragePlottR.R --normalizedCov_file ${FileNormalizedCoverage} --window_size 1000 --sample ${sample} --chr "USMA_521_v2_"

#
