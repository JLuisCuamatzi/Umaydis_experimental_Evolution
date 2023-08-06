# author: jcuamatzi
# first, clean env
rm(list = ls())
#load libraries
library_names <- c("vcfR", "dplyr", "tidyr", "tidyverse", "data.table", "ComplexHeatmap", "ggplot2", "Cairo", "cowplot", "ggthemes", "circlize")

# Load each library. Install it requires it
for (lib in library_names) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}



# read custom functions
source(file = "Functions_Filtering_SNPs.R")
## 
# change dir to path were are the VCF files

vcfDir <- paste0("vcfFiles/")

setwd(vcfDir)

# Create directories for outputs
dir.2save.SNPCalling <- paste0("../SNPCalling/")
if (dir.exists(dir.2save.SNPCalling)) {
  print("Directory already exists!")
} else {
  print (paste0("Directory ", dir.2save.SNPCalling, " does not exists!"))
  print ("Directory will be created")
  dir.create(dir.2save.SNPCalling)
}


# Directory to Save Plots
dir.2save.plots <- paste0("../SNPCalling/Plots/")
if (dir.exists(dir.2save.plots)) {
  print("Directory already exists!")
} else {
  print (paste0("Directory ", dir.2save.plots, " does not exists!"))
  print ("Directory will be created")
  dir.create(dir.2save.plots)
}

# Directory to Save Tables
dir.2save.tables <- paste0("../SNPCalling/Tables/")
if (dir.exists(dir.2save.tables)) {
  print("Directory already exists!")
} else {
  print (paste0("Directory ", dir.2save.tables, " does not exists!"))
  print ("Directory will be created")
  dir.create(dir.2save.tables)
}

## 1) get SNPs from all VCF
# create a lis with the name of the samples for the analysis
# list without colonies from line C in 35 gen with or without H2O2

# just colonies (15)
sample.list <- c("2021EE01", "2021EE02", "2021EE03", "2021EE04", "2021EE05", 
                 "2021EE06", "2021EE07", "2021EE08", "2021EE09", "2021EE10", 
                 "2021EE11", "2021EE12", "2021EE13", "2021EE14", "2021EE15", 
                 "2021EE16", "2021EE17", "2021EE18", "2021EE19", "2021EE20", 
                 "2021EE21", "2021EE22",
                 "2021EE25", "2021EE26", "2021EE27",
                 "2021EE47")

# get the SNPs using the "get.SNP" function

for (i in sample.list){
  df.name <- paste("df.vcf.", i, sep = "")
  
  vcf.file <- paste(i, "_VC_bcftools.vcf.gz", sep = "")
  
  df.tmp <- get.SNP(vcf.file = vcf.file)
  assign(df.name, df.tmp)
  rm(df.tmp, df.name, vcf.file)
}

# the output of this function are individual df for each sample with the snps
rm(i, sample.list) # remove i "index" and sample.list


# 2) Concatenate SNPs data.frames

my.vcf.df <- ls(pattern = "df.vcf.") # list objects with pattern 
my.vcf.df <- mget(my.vcf.df)         # get a list with each data frames
my.vcf.df <- bind_rows(my.vcf.df)    # concatenate all df into a single df
setDT(my.vcf.df) # transform to DT format for efficient processing with package: data.table 

rm(list = ls(pattern = "df.vcf.")) # remove individual dfs


# create background data frame 
df.vcf.SG200.BGI <- my.vcf.df[Sample == "2021EE01"] # 883 rows
# add a column with the info of sequencing technology
df.vcf.SG200.BGI$SeqTechnology <- paste0("DNBSeq")

# illumina
df.vcf.SG200.Ill <- my.vcf.df[Sample == "USMA_SG200"] # 962 rows
# add a column with the info of sequencing technology
df.vcf.SG200.Ill$SeqTechnology <- paste0("NextSeq")

## check qual and dp in both data sets
plot.Q.DP.SG200.BGI <- plot.qual.and.dp(vcf.df = df.vcf.SG200.BGI, 
                 qual.threshold = median(df.vcf.SG200.BGI$QUAL), 
                 dp.threshold = median(df.vcf.SG200.BGI$DP),
                 plot.title.dp = "SG200 sequenced with DNBSeq (All SNPs)", 
                 plot.title.q = "SG200 sequenced with DNBSeq (All SNPs)") # dashed lines indicates the median of QUAL and DP respectively

plot.Q.DP.SG200.Ill <- plot.qual.and.dp(vcf.df = df.vcf.SG200.Ill, 
                 qual.threshold = median(df.vcf.SG200.Ill$QUAL), 
                 dp.threshold = median(df.vcf.SG200.Ill$DP),
                 plot.title.dp = "SG200 sequenced with NextSeq (All SNPs)", 
                 plot.title.q = "SG200 sequenced with NextSeq (All SNPs)")

# Export Plots
ggsave(filename = paste0(dir.2save.plots, "Plot01.Qual_DP_SG200.BGI.png"), plot = plot.Q.DP.SG200.BGI, width = 18, height = 10, units = "in", dpi = 300)

ggsave(filename = paste0(dir.2save.plots, "Plot02.Qual_DP_SG200.Illumina.png"), plot = plot.Q.DP.SG200.Ill, width = 18, height = 10, units = "in", dpi = 300)

# concatenate both data frames
df.SG200.both <- bind_rows(df.vcf.SG200.BGI, 
                           df.vcf.SG200.Ill) # 1845 rows, means 1845 SNPs

## Upset plot for intersect 
# create list 
SG200.both.list <- split(df.SG200.both$SNP_ID, 
      df.SG200.both$SeqTechnology)
# make comb matrix
SG200.both.list <- make_comb_mat(SG200.both.list)
# plot and export
file2export <- paste0(dir.2save.plots, "Plot03.UpSet.SG200.png") # Create object with filename to save
CairoPNG(file = file2export, width = 10, height = 8, units = "in", dpi =300)
UpSet(m = SG200.both.list,
      set_order = c("DNBSeq", "NextSeq"),
      top_annotation = upset_top_annotation(SG200.both.list,
                                            annotation_name_rot = 90,
                                            annotation_name_side = "right",
                                            axis_param = list(side = "right"),
                                            add_numbers = T,
                                            numbers_rot = 90),
      right_annotation = upset_right_annotation(SG200.both.list, add_numbers = T))
dev.off()

# extract intersect between both sequencing technologies
SG200.Intersect <- df.SG200.both %>% filter(duplicated(SNP_ID)) # 801 rows

# identify unique SNPs for each sequencing technology
SG200.Unique.BGI <- df.vcf.SG200.BGI[ !(df.vcf.SG200.BGI$SNP_ID %in% 
                                          SG200.Intersect$SNP_ID ) , ] # 82 rows
SG200.Unique.Ill <- df.vcf.SG200.Ill[ !(df.vcf.SG200.Ill$SNP_ID %in% 
                                          SG200.Intersect$SNP_ID ) , ] # 161 rows

SG200.background <- bind_rows(SG200.Intersect, SG200.Unique.BGI, SG200.Unique.Ill) # 1044

# I can get the same snps if I remove duplicated SNPs
#nrow(df.SG200.both %>% filter(!duplicated(SNP_ID)))

# check Q and DP for uniques data sets
# bgi:
plot04 <- plot.qual.and.dp(vcf.df = SG200.Unique.BGI, 
                 plot.title.q = "Distribution of Q and DP for SNPs that were only identified in SG200 sequenced with BGI",
                 plot.title.dp = "Distribution of Q and DP for SNPs that were only identified in SG200 sequenced with BGI")
# illumina:
plot05 <- plot.qual.and.dp(vcf.df = SG200.Unique.Ill, 
                 plot.title.q = "Distribution of Q and DP for SNPs that were only identified in SG200 sequenced with Illumina",
                 plot.title.dp = "Distribution of Q and DP for SNPs that were only identified in SG200 sequenced with Illumina")

# Export plots using ggsave
ggsave(filename = paste0(dir.2save.plots, "Plot04.Qual_DP_SG200.Unique.BGI.png"), plot = plot04, width = 18, height = 10, units = "in", dpi = 300)

ggsave(filename = paste0(dir.2save.plots, "Plot05.Qual_DP_SG200.Unique.Illumina.png"), plot = plot05, width = 18, height = 10, units = "in", dpi = 300)

# Remove objects
rm(SG200.both.list, SG200.Intersect, SG200.Unique.BGI, SG200.Unique.Ill, 
   df.SG200.both, df.vcf.SG200.BGI, df.vcf.SG200.Ill )

# change 2021EE01 to USMA_SG200
# the sample sequenced with Illumina is genotyped as: USMA_SG200
# check the number of rows that match with USMA_SG200
nrow(my.vcf.df[Sample == "USMA_SG200"]) # result = 962
# the sample sequenced with BGI is genptyped as: 2021EE01
# check the number of rows that match with USMA_SG200
nrow(my.vcf.df[Sample == "2021EE01"]) # result = 883

# change 2021EE01 to USMA_SG200:
my.vcf.df$Sample <- gsub("2021EE01", "USMA_SG200", my.vcf.df$Sample) 

# check the number of rows for "USMA_SG200", the new result should be the sum of 883 + 962 = 1845
nrow(my.vcf.df[Sample == "USMA_SG200"]) # result = 1845

# remove rows that match with "USMA_SG200". I expected 21,630 rows
my.vcf.df <- my.vcf.df[Sample != "USMA_SG200"]

# remove background using custom function
# in this function we need a data frame with background SNPs and a dataframe with SNPs from the rest of samples 

filtered.vcf.df <- remove.background(background.sample.df = SG200.background, 
                                     variants.df = my.vcf.df) #  846 obs when I removed 2021EE01 (BGI) and USMA_SG200 (Illumina)
nrow(filtered.vcf.df)
# check qual and dp
# plot distribution in all SNPs without any filtering process
plot06.distribution.AllSNPs <- plot.qual.and.dp(vcf.df = my.vcf.df, 
                 plot.title.q = "All Samples - All SNPs",
                 plot.title.dp = "All Samples - All SNPs")

# plot distribution in SNPs without SG200 as background
plot07.distribution.NoSG200 <- plot.qual.and.dp(vcf.df = filtered.vcf.df, 
                 qual.threshold = median(filtered.vcf.df$QUAL), 
                 dp.threshold = median(filtered.vcf.df$DP),
                 plot.title.q = "SNPs without SG200 background",
                 plot.title.dp = "SNPs without SG200 background")

# filter by Q > 200
filtered.vcf.df.Q200 <- filtered.vcf.df %>% filter(QUAL > 200) # just 47 obs
nrow(filtered.vcf.df.Q200)
# plot distribution in SNPs after background and after filter by Q > 200
plot08.distribution.Q200 <- plot.qual.and.dp(vcf.df = filtered.vcf.df.Q200, 
                                    qual.threshold = 200, dp.threshold = 10, 
                                    plot.title.q = "SNPs without SG200 (Q > 200)", 
                 plot.title.dp = "SNPs without SG200 (Q > 200)")


# save plots
# use a if - else to avoid rewrite the files
file06 <- paste0(dir.2save.plots, "Plot06.Distribution_Q_and_DP_All_SNPs.png")
# Plot 06
if (file.exists(file06)){
  print("File already exists!")
} else {
  print("The file will be generated using ggsave")
  ggsave(filename = file06,
         plot = plot06.distribution.AllSNPs, width = 16, height = 8, dpi = 300, units = "in")
}

# Plot 07
file07 <- paste0(dir.2save.plots, "Plot07.Distribution_Q_and_DP_SNPs_wo_background.png")

if (file.exists(file07)){
  print("File already exists!")
} else {
  print("The file will be generated using ggsave")
  ggsave(filename = file07,
         plot = plot07.distribution.NoSG200, width = 16, height = 8, dpi = 300, units = "in")
}

# Plot 08
file08 <- paste0(dir.2save.plots, "Plot08.Distribution_Q_and_DP_SNPs_wo_background.Q200.png")
if (file.exists(file08)){
  print("File already exists!")
} else {
  print("The file will be generated using ggsave")
  ggsave(filename = file08,
         plot = plot08.distribution.Q200, width = 16, height = 8, dpi = 300, units = "in")
}

# remove objects that will no longer be used

rm(plot06.distribution.AllSNPs, plot07.distribution.NoSG200, plot08.distribution.Q200)


### work with SNPs in Q > 200
df.info <- fread("../../SampleInfo.csv")
# add sample info to filtered vcf q >200
filtered.vcf.df.Q200 <- filtered.vcf.df.Q200 %>% 
  left_join(select(df.info, ID, Name, Line, OrderName), by = c("Sample" = "ID"))
# re-order data frame by OrderName
filtered.vcf.df.Q200 <- filtered.vcf.df.Q200 %>% arrange(OrderName)

filtered.vcf.df.Q200$SNP_ID_Sample <- paste0(filtered.vcf.df.Q200$SNP_ID, "_", 
                                             filtered.vcf.df.Q200$Sample)

# export this table: filtered.vcf.df.Q200
write.table(x = filtered.vcf.df.Q200, file = paste0(dir.2save.tables,"Filtered.Variants.Q200.csv"), sep = ",", row.names = F, quote = F)

#  count SNP frequency
SNP.Q200.Frequency <- filtered.vcf.df.Q200 %>% group_by(SNP_ID) %>% summarise(SNP.Frequency = n()) %>% ungroup() %>% arrange(-SNP.Frequency)

shared.SNPs.Q200 <- filtered.vcf.df.Q200 %>% 
  select(!c("QUAL", "DP", "Sample", "Name", "Line", "SNP_ID_Sample")) %>% 
  #group_by(SNP_ID) %>% 
  pivot_wider(names_from = OrderName, 
              values_from = GT)
# change NA by 0
shared.SNPs.Q200[is.na(shared.SNPs.Q200)] <- 0

# add SNP frequency to df with all snps (shared.SNPs.Q200)
shared.SNPs.Q200 <- shared.SNPs.Q200  %>% 
  left_join(select(SNP.Q200.Frequency, SNP_ID, SNP.Frequency),
            by = c("SNP_ID" = "SNP_ID"))
# move column '33' between '5' and '6'
shared.SNPs.Q200 <- shared.SNPs.Q200[, c(1:5, 23, 6:22)]
# re order df based on SNP frequency 
shared.SNPs.Q200 <- shared.SNPs.Q200 %>% arrange(-SNP.Frequency)

## Export Tables of Shared SNPs and their Frequency in the Samples (shared.SNPs.Q200)
write.table(x = shared.SNPs.Q200, file = paste0(dir.2save.tables, "Shared.SNP.NoSG200.Q200.csv"), sep = ",", quote = F, row.names = F)

## Create BED File to Verify Allele Frequence
# BED file must have three columns. 1.- Chr, 2.- (SNP Position - 1), 3.- SNP Position. BED File must not have Header
snps.bed.file <- data.frame("C1" = shared.SNPs.Q200$CHROM, "C2" = (shared.SNPs.Q200$POS - 1), "C3" = shared.SNPs.Q200$POS)
# export
write.table(x = snps.bed.file, file = paste0(dir.2save.tables, "SNPsToCheck.bed"), sep = "\t", quote = F, row.names = F, col.names = F )

# Count SNP per Sample in each Filtering Step
SNP.per.Sample.0 <- my.vcf.df %>% group_by(Sample) %>% summarise(Num.SNP_All = n())
SNP.per.Sample.1 <- filtered.vcf.df %>% group_by(Sample) %>% summarise(Num.SNP_NoSG200 = n())
SNP.per.Sample.2 <- filtered.vcf.df.Q200 %>% group_by(Sample) %>% summarise(Num.SNP_NoSG200.Q200 = n())

SNP.Per.Sample <- SNP.per.Sample.0 %>% 
  left_join(select(SNP.per.Sample.1, Sample, Num.SNP_NoSG200), by = c("Sample" = "Sample"))
SNP.Per.Sample <- SNP.Per.Sample %>% 
  left_join(select(SNP.per.Sample.2, Sample, Num.SNP_NoSG200.Q200), by = c("Sample" = "Sample"))
# Convert NA to 0
SNP.Per.Sample[is.na(SNP.Per.Sample)] <- 0
# Add Info
SNP.Per.Sample <- SNP.Per.Sample %>% 
  left_join(select(df.info, ID, Name, Line), by = c("Sample" = "ID"))
# Export table
write.table(x = SNP.Per.Sample, file = paste0(dir.2save.tables,"SNP_Number_Per_Sample.Q200.csv"), sep = ",", row.names = F, quote = F)


## Create Matrix of Shared SNPs
# create list with SNP ID and Name (sample name)
SNP.List.Q200 <- split(filtered.vcf.df.Q200$SNP_ID, 
                       factor(filtered.vcf.df.Q200$OrderName))
# create matrix raw and percentage
Mtrx.Q200.Raw <- intersect.snp.raw(list_snp = SNP.List.Q200)
Mtrx.Q200.Per <- intersect.snp.perc.NoSym(list_snp = SNP.List.Q200) # no symmetric matrix

#export matrix as csv
# matrix raw
write.table(Mtrx.Q200.Raw, file = paste0(dir.2save.tables,"Matriz.Q200.Raw.csv"), sep = ",", quote = F)
write.table(Mtrx.Q200.Per, file = paste0(dir.2save.tables,"Matriz.Q200.Percentage.csv"), sep = ",", quote = F)

# function to col the heat map
col_fun <- colorRamp2(c(0, max(Mtrx.Q200.Raw)/2, max(Mtrx.Q200.Raw)), 
                      c("yellowgreen", "seagreen4", "navy")) # colors for raw matrix

# create groups to split the matrix (col and row)
names(SNP.List.Q200)
row_split <- c(rep("A", length(grep("LA", names(SNP.List.Q200)))),
               rep("B", length(grep("LB", names(SNP.List.Q200)))),
               rep("C", length(grep("LC", names(SNP.List.Q200)))),
               rep("D", length(grep("LD", names(SNP.List.Q200)))),
               rep("E", length(grep("LE", names(SNP.List.Q200)))),
               rep("F", length(grep("LF", names(SNP.List.Q200))))
               )

row_split <- factor(row_split, 
                    levels = c("A", "B", "C", "D", "E", "F"))

## export matrix using cairo package
file2export <- paste0(dir.2save.plots, "Matrix.Q200.RAW.svg") # Create object with filename to save
CairoSVG(file = file2export, width = 20, height = 15)
Heatmap(Mtrx.Q200.Raw,
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
          grid.text(Mtrx.Q200.Raw[i, j], x, y)
        })
dev.off()

# percentage
col_fun <- colorRamp2(c(0, 50, 100), 
                      c("yellowgreen", "seagreen4", "navy"))
Mtrx.Q200.Per <- (Mtrx.Q200.Per*100) # multiply per 100 

Mtrx.Q200.Per <- as.data.frame(Mtrx.Q200.Per) %>% 
  mutate(across(where(is.numeric), ~ round(., 0))) # round the percentage to integer 

Mtrx.Q200.Per <- as.matrix(Mtrx.Q200.Per) # turn into matrix 

# export
file2export <- paste0(dir.2save.plots, "Matrix.Q200.PER.svg") # Create object with filename to save
CairoSVG(file = file2export, width = 20, height = 15)
Heatmap(Mtrx.Q200.Per,
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
          grid.text(Mtrx.Q200.Per[i, j], x, y)
        })
dev.off()

# script ends
rm(list = ls() ) # clean env
# This Script Ends in this Point

# The next step is verify the Allele Frequency (of the SNPs identified with Q > 200) in the pools and also, verify its coverage in the colonies

