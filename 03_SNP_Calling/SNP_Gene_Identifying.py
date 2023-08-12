

"""
@ author: j cuamatzi

script to identify if a SNP is in a gene (ORF) or it is intergenic

@ usage:

cd ~/03_SNP_Calling/

python3 SNP_Gene_Identifying.py -s SNPCalling/Tables/Shared.SNP.NoSG200.Q200.AF90.csv -g ../USMA_521_GeneProteines_DB.csv -o Umaydis_EE_Annotated_SNPs.csv


"""


import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description="Detected if the SNP is in the ORF of a gene or is intergenic")
parser.add_argument("-s", "--snps_file", help="CSV file containing the list of SNPs", required=True)
parser.add_argument("-g", "--genes_database", help="CSV file containing the coordinates of the umaydis genes", required=True)
parser.add_argument("-o", "--output_file", help="Name for the output file (it should contain the csv extension)", required=True)
args = parser.parse_args() 

# Read file with SNPs
snps_data = pd.read_csv(args.snps_file)
new_column_names = {'CHROM': 'Chromosome', 'POS': 'Position', 'REF': 'Reference', 'ALT':'Alternative'}
snps_data.rename(columns=new_column_names, inplace=True)

snps_data_1 = snps_data.iloc[:, 0:2]


# Read U. maydis gene data base
usma_genes = pd.read_csv(args.genes_database)
usma_genes_1 = pd.concat([usma_genes.iloc[:, 0:5], usma_genes.iloc[:, 8:9]], axis=1)

# Create an empty result list
result_list = []

# Iterate through snps_data_1
for index, row in snps_data_1.iterrows():
    matching_gene = usma_genes_1[(usma_genes_1["Chromosome"] == row["Chromosome"]) &
                                 (usma_genes_1["Start_position"] <= row["Position"]) &
                                 (usma_genes_1["End_position"] >= row["Position"])]
    
    if not matching_gene.empty:
        gene_info = matching_gene.iloc[0]
        result_list.append({"Chromosome": row["Chromosome"],
                            "Position": row["Position"],
                            "Gene": gene_info["UMAG_ID"],
                            "Start_position": gene_info["Start_position"],
                            "End_position": gene_info["End_position"],
                            "ProteinName": gene_info["ProteinNames"]})
    else:
        result_list.append({"Chromosome": row["Chromosome"],
                            "Position": row["Position"],
                            "Gene": "Intergenic",
                            "Start_position": "NA",
                            "End_position": "NA",
                            "ProteinName": "NA"})

# Create the result DataFrame from the list
result = pd.DataFrame(result_list)

# Insert the new columns after 'SNP_ID'
result = pd.concat([snps_data.iloc[:, 1:5], result.iloc[:,0:6], snps_data.iloc[:, 5:23]], axis = 1)

# Export "result" as csv
result.to_csv(args.output_file, index=False)  
