# author: jcuamatzi
# first, clean env
rm(list = ls())

# Load libraries
# Libraries
library_names <- c("vcfR", "dplyr", "tidyr", "tidyverse", "data.table", "ComplexHeatmap", "ggplot2", "Cairo", "cowplot", "ggthemes", "circlize", "scales")

# Load each library. Install it requires it
for (lib in library_names) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

# Create directoris to export plots and tables

  # Directory to Save Plots
dir.2save.plots <- paste0( "SNPCalling/Plots/")
if (dir.exists(dir.2save.plots)) {
  print (paste0("Directory ", dir.2save.plots, " already exists!"))
} else {
  print (paste0("Directory ", dir.2save.plots, " does not exists!"))
  print (paste0("Directory ", dir.2save.plots, " will be created!"))
  dir.create(dir.2save.plots)
}

# Directory to Save Tables
dir.2save.tables <- paste0("SNPCalling/Tables/")
if (dir.exists(dir.2save.tables)) {
  print (paste0("Directory ", dir.2save.tables, " already exists!"))
} else {
  print (paste0("Directory ", dir.2save.tables, " does not exists!"))
  print (paste0("Directory ", dir.2save.tables, " will be created!"))
  dir.create(dir.2save.tables)
}

# Read custom functions to process SNPs
source(file = "Functions_Filtering_SNPs.R")

# Read sample sheet and files with more metadata from samples
df.samples <- fread("../USMA_EE_Colonies_SampleSheet.csv") # sample name
df.info <- fread("../SampleInfo.csv") # metadata about generations, treatment group, etc

list.colonies <- df.samples$SampleID

# change current directory to directory that contains files: {sampleID}_SNPs_AlleleDepth.txt
setwd("SNPCalling/Tables/Check_SNPs_AF/")
# A) Read .txt files with depth allele coverage
# use the next for loop to read the files and order the columns into a data frame
for (i in list.colonies) {
  df.name <- paste0("df.colony.", i)
  
  #file2read <- paste0(i, "_SNPs_in_LineC.txt") # snps indentified only in cols at 200 line C
  file2read <- paste0(i, "_SNPs_AlleleDepth.txt")
  
  df.tmp <- fread(file2read)
  
  names(df.tmp) <- c("Chr", "Position", "Ref", "Alt", "Depth")  # rename cols
  df.tmp$Sample <- paste0(i)    # here, I created a column to paste the sample
  
  
  df.tmp <- transform.table.snp.freq(df.snp.freq = df.tmp)
  
  assign(df.name, df.tmp)
  rm(df.tmp, i, file2read, df.name)
  
}

# B) Concatenate the data frame generated in the previous for-loop into a single df
df.Colony.Frequency <- ls(pattern = "df.colony.") # list objects with pattern 
df.Colony.Frequency <- mget(df.Colony.Frequency)         # get a list with each data frames
df.Colony.Frequency <- bind_rows(df.Colony.Frequency)    # concatenate all df into a single df

rm(list = ls(pattern = "df.colony.2021")) # remove individual dfs

setwd("../../../") # go back to the directory "~/03_SNP_Calling/"

# add line information from df.samples to 'df.Colony.Frequency'
df.Colony.Frequency <- df.Colony.Frequency %>% left_join(select(df.info, ID, Name, Line), by = c("Sample" = "ID"))

#read SNPs info (SNPs identified after Q200)
df.SNP.Q200 <- fread("SNPCalling/Tables/Shared.SNP.NoSG200.Q200.csv")

df.SNP.Q200.tmp <- df.SNP.Q200

df.SNP.Q200.tmp <- df.SNP.Q200.tmp %>% select( CHROM, POS, ALT, SNP_ID) %>% mutate(SNP_ID_1 = paste0(CHROM, "-", POS))
names(df.SNP.Q200.tmp) <- c("Chr", "Pos", "Alt_Detected", "SNP_ID", "SNP_ID_1")

# add SNP Information (SNP_ID, SNP_ID_1 and Alt_Detected)
df.Colony.Frequency <- df.Colony.Frequency %>% 
  left_join(select(df.SNP.Q200.tmp, SNP_ID_1, SNP_ID, Alt_Detected),
                                  by = c("SNP_ID_1" = "SNP_ID_1"))
# create a string that joins SNP_ID and Sample
df.Colony.Frequency$SNP_ID_Sample <- paste0(df.Colony.Frequency$SNP_ID, "_",
                                            df.Colony.Frequency$Sample)

# Modify col 'NucleotideType'
df.Colony.Frequency <- df.Colony.Frequency %>% group_by(SNP_ID) %>% 
  mutate(NucleotideType = if_else(Nucleotide == Alt_Detected, 
                                  paste0("Alt (", Alt_Detected, ")"),
                                  NucleotideType),
         NucleotideType = if_else(Nucleotide != Alt_Detected & 
                                    Nucleotide != Ref,
                                  "Other", NucleotideType)) %>% 
  ungroup()

# remove SNP_ID_1
df.Colony.Frequency <- df.Colony.Frequency %>% select(!SNP_ID_1)

# create df to concatenate with filtered.vcf.df.Q200 # in 'df.Colony.Frequency$NucleotideType' remove Reference and others
df.2concatenate <- df.Colony.Frequency %>% 
  filter(!NucleotideType %in% c("Reference", "Other") ) %>% 
  select(NucleotideType, PerFreq, SNP_ID_Sample, Sample)

# read filtered df q200 from step 1  
filtered.vcf.df.Q200 <- fread(  "SNPCalling/Tables/Filtered.Variants.Q200.csv")

filtered.vcf.df.Q200 <- filtered.vcf.df.Q200 %>% 
  left_join(select(df.2concatenate, NucleotideType, PerFreq, SNP_ID_Sample), 
            by = c("SNP_ID_Sample" = "SNP_ID_Sample"))  
  
# filter by threshold of PerFreq (Allele frequency in percentage) >= 0.9 (meaning 90%)
# 
filtered.vcf.df.Q200.AF90 <- filtered.vcf.df.Q200 %>% 
  filter(PerFreq >= 0.9) # reduce to 45 obs from 47
# count SNP frequency
SNP.Num.Q200.AF90 <- filtered.vcf.df.Q200.AF90 %>% 
  group_by(SNP_ID) %>% 
  summarise(SNP.Frequency.Q200.AF90 = n())
  
# remove SNPs in SG200 that have an Allele Frequency > 40%
# here: Extract Allele depth in SG200 (both samples)
SNPs.in.SG200 <- df.Colony.Frequency %>% 
  filter(Sample %in% c("2021EE01", "2021EE47"))
# filter by 0.4 (meaning 40%):
SNPs.in.SG200 <- SNPs.in.SG200 %>% filter(PerFreq > 0.4)
# remove rows with 'Reference' or 'Other' in 'NucleotideType'
SNPs.in.SG200 <- SNPs.in.SG200 %>% filter(!NucleotideType %in% c("Reference", "Other"))

# remove SNPs in SG200 that have an Allele Frequency > 40%
filtered.vcf.df.Q200.AF90 <- filtered.vcf.df.Q200.AF90[ c(!filtered.vcf.df.Q200.AF90$SNP_ID 
                                                          %in% 
                                                            SNPs.in.SG200$SNP_ID), ]  # from 45 to 43
# count SNP freq.
SNP.Num.Q200.AF90.NoSG200AF40 <- filtered.vcf.df.Q200.AF90 %>% 
  group_by(SNP_ID) %>% 
  summarise(SNP.Frequency.Q200.AF90.NS = n())

# prepare table for shared SNPs
# select just the next columns
shared.SNPs.Q200.AF90 <- filtered.vcf.df.Q200.AF90 %>% 
  select(!c("QUAL", "DP", "Sample", "Name", "Line", "SNP_ID_Sample", "NucleotideType", "PerFreq")) %>% 
  #group_by(SNP_ID) %>% 
  pivot_wider(names_from = OrderName, 
              values_from = GT)
# change NA by 0
shared.SNPs.Q200.AF90[is.na(shared.SNPs.Q200.AF90)] <- 0
# count SNP frequency
SNP.Freq.Q200.AF90 <- filtered.vcf.df.Q200.AF90 %>% 
  group_by(SNP_ID) %>% count(name = "SNP.Frequency.Q200.AF90")
# add SNP frequency to df with all snps (shared.SNPs.Q200)
shared.SNPs.Q200.AF90 <- shared.SNPs.Q200.AF90  %>% 
  left_join(select(SNP.Freq.Q200.AF90, SNP_ID, SNP.Frequency.Q200.AF90),
            by = c("SNP_ID" = "SNP_ID"))
  
  # move column '33' between '5' and '6'
shared.SNPs.Q200.AF90 <- shared.SNPs.Q200.AF90[, c(1:5, 23, 6:22)]
# re order df based on SNP frequency 
shared.SNPs.Q200.AF90 <- shared.SNPs.Q200.AF90 %>% arrange(-SNP.Frequency.Q200.AF90)

# export table
write.table(x = shared.SNPs.Q200.AF90, file = paste0(dir.2save.tables, "Shared.SNP.NoSG200.Q200.AF90.csv"), sep = ",", quote = F, row.names = F)

# count SNP per sample after this filtering
SNP.per.Sample.4 <- filtered.vcf.df.Q200.AF90 %>% group_by(Sample) %>% summarise(Num.SNP_NoSG200.Q200.AF90 = n())

SNP.per.Sample <- fread("SNPCalling/Tables/SNP_Number_Per_Sample.Q200.csv")
# Concatenate SNP.per.Sample.4 in SNP.per.Sample
SNP.per.Sample <- SNP.per.Sample %>% left_join(select(SNP.per.Sample.4, Sample, Num.SNP_NoSG200.Q200.AF90), by = c("Sample" = "Sample"))
# change NA by 0
SNP.per.Sample[is.na(SNP.per.Sample)] <- 0
SNP.per.Sample <- SNP.per.Sample[,c(1:4, 7, 5:6)]
# Export table
write.table(x = SNP.per.Sample, file = paste0(dir.2save.tables,"SNP_Number_Per_Sample.Q200.AF90.csv"), sep = ",", row.names = F, quote = F)

## Create Matrix of Shared SNPs
# create list with SNP ID and Name (sample name)
SNP.List.Q200.AF90 <- split(filtered.vcf.df.Q200.AF90$SNP_ID, 
                       factor(filtered.vcf.df.Q200.AF90$OrderName))
# create matrix raw and percentage
Mtrx.Q200.AF90.Raw <- intersect.snp.raw(list_snp = SNP.List.Q200.AF90)
Mtrx.Q200.AF90.Per <- intersect.snp.perc.NoSym(list_snp = SNP.List.Q200.AF90) # no symmetric matrix

#export matrix as csv
# matrix raw
write.table(Mtrx.Q200.AF90.Raw, file = paste0(dir.2save.tables,"Matriz.Q200.AF90.Raw.csv"), sep = ",", quote = F)
write.table(Mtrx.Q200.AF90.Per, file = paste0(dir.2save.tables,"Matriz.Q200.AF90.Percentage.csv"), sep = ",", quote = F)

# function to col the heat map
col_fun <- colorRamp2(c(0, max(Mtrx.Q200.AF90.Raw)/2, max(Mtrx.Q200.AF90.Raw)), 
                      c("yellowgreen", "seagreen4", "navy")) # colors for raw matrix

# create groups to split the matrix (col and row)
names(SNP.List.Q200.AF90)
row_split <- c(rep("A", length(grep("LA", names(SNP.List.Q200.AF90)))),
               rep("B", length(grep("LB", names(SNP.List.Q200.AF90)))),
               rep("C", length(grep("LC", names(SNP.List.Q200.AF90)))),
               rep("D", length(grep("LD", names(SNP.List.Q200.AF90)))),
               rep("E", length(grep("LE", names(SNP.List.Q200.AF90)))),
               rep("F", length(grep("LF", names(SNP.List.Q200.AF90))))
               )

row_split <- factor(row_split, 
                    levels = c("A", "B", "C", "D", "E", "F"))

## export matrix using cairo package
file2export <- paste0(dir.2save.plots, "Matrix.Q200.AF90.RAW.svg") # Create object with filename to save
CairoSVG(file = file2export, width = 20, height = 15)
Heatmap(Mtrx.Q200.AF90.Raw,
        name = "Shared SNP\n(Raw Numbers)", 
        cluster_rows = F, 
        cluster_columns = F,
        column_split = row_split,
        row_split = row_split,
        border = T,
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        col = col_fun,
        cell_fun=function(j, i, x, y, w, h, col){
          grid.text(Mtrx.Q200.AF90.Raw[i, j], x, y)
        })
dev.off()

# percentage
col_fun <- colorRamp2(c(0, 50, 100), 
                      c("yellowgreen", "seagreen4", "navy"))
Mtrx.Q200.AF90.Per <- (Mtrx.Q200.AF90.Per*100) # multiply per 100 

Mtrx.Q200.AF90.Per <- as.data.frame(Mtrx.Q200.AF90.Per) %>% 
  mutate(across(where(is.numeric), ~ round(., 0))) # round the percentage to integer 

Mtrx.Q200.AF90.Per <- as.matrix(Mtrx.Q200.AF90.Per) # turn into matrix 

# export
file2export <- paste0(dir.2save.plots, "Matrix.Q200.AF90.PER.svg") # Create object with filename to save
CairoSVG(file = file2export, width = 20, height = 15)
Heatmap(Mtrx.Q200.AF90.Per,
        name = "Shared SNP\n(%)", 
        cluster_rows = F, 
        cluster_columns = F,
        column_split = row_split,
        row_split = row_split,
        border = T,
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        col = col_fun,
        cell_fun=function(j, i, x, y, w, h, col){
          grid.text(Mtrx.Q200.AF90.Per[i, j], x, y)
        })
dev.off()

# script ends
rm(list = ls() ) # clean env
# This Script Ends in this Point

# The next step is verify the Allele Frequency (of the SNPs identified with Q > 200) in the pools and also, verify its coverage in the colonies

