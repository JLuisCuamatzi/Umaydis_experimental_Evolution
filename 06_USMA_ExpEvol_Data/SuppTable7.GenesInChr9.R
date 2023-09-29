# author: jcuamatzi
# Script to identify the genes located in the first 160 kb of chr 9
# supp table 7

## Load libraries
libraries <- c("data.table", "dplyr", "writexl", "openxlsx")


# Check: if libraries is does not install, install, then just load it
for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

rm(list = ls())

# Set as working directory, the directory in which this script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read file with USMA 521 v2 genes
df <- fread("../USMA_521_GeneProteines_DB.csv", fill = T)

# Extract Chr9
df.Chr9 <- df %>% filter(Chromosome == "USMA_521_v2_9")

# Check the first 160 kb
df.Chr9.160kb <- df.Chr9 %>% filter(End_position < 160000) %>% 
  select(2,3,4,5, 9)


# check if the genes are within the amplification, use as threshold the putative breakpoint identified by CNVnator (149100)

df.Chr9.160kb <- df.Chr9.160kb %>% mutate(
  InAmplification = if_else(End_position < 149100, "Yes", "No")
)


# Change headers to export
names(df.Chr9.160kb) <- c("Chr", "Start", "End", "GeneID", "ProteinName", "InAmplification")

# remove USMA_521_v2_
df.Chr9.160kb$Chr <- gsub("USMA_521_v2_", "", df.Chr9.160kb$Chr)

# Export as word document

# doc <- read_docx()
# 
# table2save <- flextable(df.Chr9.160kb)
# 
# doc <- body_add_flextable(doc, value = table2save)
# 
# 
# print(doc, target = "SuppTable7.GenesIn160kbChr9.docx")
# file.copy("SuppTable7.GenesIn160kbChr9.docx", getwd())  


## Update: export as Excel
# df.Chr9.160kb
excel_file <- "Supp.Table7.GenesIn160kbChr9.xlsx"

write.xlsx(x = df.Chr9.160kb, file = excel_file)

rm(list = ls())








