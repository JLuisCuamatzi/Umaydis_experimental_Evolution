
# script for figure 4 in the main manuscript
# author: jcuamatzi

# Libraries
library_names <- c("data.table", "ggthemes", "ggplot2", "tidyverse", "dplyr", "scales", "cowplot", "grid", "optparse")

# Load each library. Install it requires it
for (lib in library_names) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

rm(lib, library_names)

# Custom functions:
transform.table.snp.freq <- function(df.snp.freq){
  df_split <- data.frame() # create empty df
  
  # work with original df
  df.snp.freq$Nucleotide <- paste0(df.snp.freq$Ref, ",", df.snp.freq$Alt) # first nt is reference
  df.snp.freq <- df.snp.freq[,c(1:3,7,5:6)]     # reorder df
  
  # for loop to pivot the data
  for (i in seq_len(nrow(df.snp.freq))) {
    nt_list <- unlist(strsplit(df.snp.freq$Nucleotide[i], ","))
    depth_list <- unlist(strsplit(df.snp.freq$Depth[i], ","))
    
    num_rows <- max(length(nt_list), length(depth_list))
    
    nt_list <- rep(nt_list, length.out = num_rows)
    
    depth_list <- rep(depth_list, length.out = num_rows)
    
    df_split <- rbind(df_split, data.frame(
      Chr = rep(df.snp.freq$Chr[i], num_rows),
      Position = rep(df.snp.freq$Position[i], num_rows),
      Ref = rep(df.snp.freq$Ref[i], num_rows),
      Nucleotide = nt_list,
      Depth = as.numeric(depth_list),
      Sample = rep(df.snp.freq$Sample[i], num_rows) ))
  }
  
  # end for loop
  df_split <- df_split %>% 
    mutate(NucleotideType = if_else(Nucleotide == Ref, "Reference", "Alternative"))
  
  df_split$SNP_ID_1 <- paste0(df_split$Chr, "-", df_split$Position)
  
  df_split <- df_split %>% group_by(SNP_ID_1) %>% 
    mutate(PerFreq = Depth/sum(Depth)) %>% 
    ungroup()
   # return pivoted df
  return(df_split)
}


## function to adjust the coverage
adjust.coverage <- function(dataset.coverage, threshold.cov){
  
  threshold.cov <- as.numeric(threshold.cov)
  
  
  
  # estimate normalized coverage (coverage in each window/global median coverage)
  dataset.coverage <- dataset.coverage %>% 
    group_by(Sample) %>% 
    # Estimate a global normalized coverage for each sample
    mutate(Global.Norm.Cov = median(normalized.coverage)) %>% 
    
    
    mutate(Adj.Norm.Cov = if_else(normalized.coverage <= (Global.Norm.Cov * threshold.cov), 
                                  0.1,
                                  normalized.coverage)) %>% 
    ungroup() %>% setDT()
  
  # select a specific chromosome
  #df.cov.chr <- dataset.coverage %>% filter(chr == target.chr)
  
  
  return(dataset.coverage)  

}

# Check Duplication at Chr 9

setwd("../02_MappingAndCoveragePlotting/normalizedCoverageTables/")

## Read Files with Normalized Coverage 
sample.list <- c("2021EE30", "2021EE31", "2021EE32", "2021EE33",
                 "2021EE34", "2021EE35", "2021EE36", "2021EE36",
                 "2021EE01") # Define list with samples

# For loop to read the files
for (sample in sample.list) {
  file.2read <- paste0(sample, "_NormalizedCoverage.txt")
  df.name <- paste0("df.", sample)
  
  df.tmp <- fread(file.2read)
  
  df.tmp$Sample <- sample
  
  assign(df.name, df.tmp)
  
}

rm(df.name, df.tmp)

cov.df <- bind_rows(mget(ls(pattern = "df.")))

setDT(cov.df)

rm(list = ls(pattern = "df."))


# Adjust the normalized coverage
df.cov.adj <- adjust.coverage(dataset.coverage = cov.df, threshold.cov = 0.25 )

# Pivot the table (from long to wide)
df.cov.adj <- df.cov.adj %>% 
  select(chr, window.start, window.end, Sample, Adj.Norm.Cov) %>% 
  pivot_wider(names_from = Sample, values_from = Adj.Norm.Cov)
   
# Estimate the Ratio of each sample / SG200
df.cov.adj.ratio <- df.cov.adj %>% 
  mutate(across(c(`2021EE01`:`2021EE36`), 
                .fns = ~./`2021EE01`, 
                .names = "{.col}_Ratio"))

# Pivot from Wide to Long
df.cov.adj.ratio <- df.cov.adj.ratio %>% 
  select(chr, window.start, window.end, `2021EE01_Ratio`:`2021EE36_Ratio`) %>% 
  pivot_longer(cols = `2021EE01_Ratio`:`2021EE36_Ratio`, 
               names_to = "Sample", values_to = "Ratio")

# Remote "_Ratio" in the column = 'Sample'
df.cov.adj.ratio$Sample <- gsub("_Ratio", "", df.cov.adj.ratio$Sample)

# Extract chromosome 9 = "USMA_521_v2_9"
setDT(df.cov.adj.ratio)
df.cov.9.pools <- df.cov.adj.ratio[chr == "USMA_521_v2_9"]

# Create labels for the x.axis
df.cov.9.pools <- df.cov.9.pools %>% mutate(x.axis = 
                       case_when(
                         startsWith("2021EE01_SG200", df.cov.9.pools$Sample) ~ "SG200",
                         startsWith("2021EE30", df.cov.9.pools$Sample) ~ "20",
                         startsWith("2021EE31", df.cov.9.pools$Sample) ~ "30",
                         startsWith("2021EE32", df.cov.9.pools$Sample) ~ "50",
                         startsWith("2021EE33", df.cov.9.pools$Sample) ~ "70",
                         startsWith("2021EE34", df.cov.9.pools$Sample) ~ "100",
                         startsWith("2021EE35", df.cov.9.pools$Sample) ~ "140",
                         startsWith("2021EE36", df.cov.9.pools$Sample) ~ "200" ))

# Indicate levels for that labels
df.cov.9.pools$x.axis <- factor(df.cov.9.pools$x.axis, 
                                levels = c(
                                  "SG200", "20", 
                                  "30",  "50", 
                                  "70",  "100", 
                                  "140", "200"))

# Estimate the Log2Ratio
df.cov.9.pools$log2Ratio <- log2(df.cov.9.pools$Ratio)

# Plot
Figure.Chr9.Pools <- df.cov.9.pools %>% 
  filter(!Sample %in% c("2021EE01") )%>% 
  ggplot(aes(x = (window.end/1000), y = log2Ratio)) +
  
  #geom_area(color = "black", fill = "gray90") +
  facet_grid(rows = vars(x.axis), scales = "free") +
  # add horizontal lines:
  geom_hline(yintercept = -2, linetype = "dotted", color = "darkgrey", linewidth = 0.3)+ 
  geom_hline(yintercept = -1, linetype = "dotted", color = "grey", linewidth = 0.3)+ 
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6)+ 
  geom_hline(yintercept = 1, linetype = "dotted", color = "darkgrey", linewidth = 0.3)+ 
  geom_hline(yintercept = 2, linetype = "dotted", color = "grey", linewidth = 0.3)+ 
  geom_hline(yintercept = 3, linetype = "dotted", color = "grey", linewidth = 0.3)+ 
  geom_line(linewidth = 0.4, color = "brown") +
  theme_classic() +
  scale_y_continuous(limits = c(-1, 3), breaks = c( 0, 1, 2, 3),
                     sec.axis = sec_axis(~./1+1, name = "Generation"))+
  scale_x_continuous(limits = c(0, max(df.cov.9.pools$window.end)/1000),
                     breaks = c(0, 150, 300, 450, 600, 733)) +
  
  labs(x = "\nChromosome 9 (kb)", y = expression("log"["2"]*" (Ratio)"),
       title =  "Coverage of Chromosome 9 in the Line C \n") +
  theme(
    axis.line.x = element_blank(),
    #axis.line.y = element_blank(),
    # titles
    #plot.title = element_text(size = 12, color = "black", hjust = 0.5),
    plot.title = element_text(hjust = 0.5, size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    # axis text
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    #
    # Second y-Axis
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    # strip
    #strip.background = element_blank(),
    strip.text = element_text(size = 10, color = "black"),
    
    panel.spacing = unit(0.5, "cm", data = NULL) )
  
# Edit Figure 4 B
Figure.Chr9.Pools <- ggplotGrob(Figure.Chr9.Pools)

lg <- linesGrob(x=unit(c(0,0),"npc"), y=unit(c(0,1),"npc"), 
                gp=gpar(col="black", lwd=4))

for (k in grep("strip-r", Figure.Chr9.Pools$layout$name)) {
  Figure.Chr9.Pools$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}


grid.draw(Figure.Chr9.Pools)



###
###
###
##### Analysis of SNP in UMAG_05545 (Chr18:330272)

setwd("../../04_Pools_analysis/")   # move to directory

# Create list of samples
list.pools <- c("2021EE30", "2021EE31", "2021EE32", "2021EE33", "2021EE34", "2021EE35", "2021EE36")

# Read files with allele frequency at that position
for (i in list.pools) {
  df.name <- paste0("df.pools.", i)
  
  #file2read <- paste0(i, "_SNPs_in_LineC.txt") # snps indentified only in cols at 200 line C
  file2read <- paste0(i, "_SNP_UMAG_05545.txt")
  
  df.tmp <- fread(file2read, sep = "\t")
  
  names(df.tmp) <- c("Chr", "Position", "Ref", "Alt", "Depth")
  df.tmp$Sample <- paste0(i)
  
  
  df.tmp <- transform.table.snp.freq(df.snp.freq = df.tmp)
  
  assign(df.name, df.tmp)
  rm(df.tmp, i, file2read, df.name)
  
}

# Concatenate the dfs into a single df
df.pools.Frequency <- bind_rows(mget(ls(pattern = "df.pools."))) # list objects with pattern and concatenate them
setDT(df.pools.Frequency)
rm(list = ls(pattern = "df.pools.2021"))

# Read the final SNPs to check the alternative detected at chr18:330272
df.SNP.ID <- fread("../03_SNP_Calling/Annotation_of_SNPs.csv") # snps (Q > 200 & AF > 90)
df.SNP.ID$SNP_ID_1 <- paste0(df.SNP.ID$Chromosome, "-", df.SNP.ID$Position)

# add whole SNP Information 
df.SNP.Target <- df.pools.Frequency %>% left_join(select(df.SNP.ID, SNP_ID_1, Alternative),by = c("SNP_ID_1" = "SNP_ID_1"))

# Paste the alternative
df.SNP.Target <- df.SNP.Target %>% mutate(NucleotideType = if_else((Nucleotide == Alternative), 
                                                     paste0("Alternative (", Alternative, ")") ,
                                                     NucleotideType))
# Paste the others Nt
df.SNP.Target <- df.SNP.Target  %>% 
  mutate(NucleotideType = if_else(Nucleotide != Alternative & Nucleotide != Ref, "Other (G/T)",
                                  NucleotideType))
  
# Paste the reference
df.SNP.Target <- df.SNP.Target %>% mutate(NucleotideType = if_else((Nucleotide == Ref), 
                                                                   paste0("Reference (", Ref, ")") ,
                                                                   NucleotideType))

# Set levels
df.SNP.Target$NucleotideType <- factor(df.SNP.Target$NucleotideType, 
                                       levels = c("Reference (C)", "Alternative (A)", "Other (G/T)"))

## Plotting the Frequency
Figure.SNPs.Pools<- df.SNP.Target %>% 
  ggplot(aes(x = Sample, y = PerFreq, 
             color = NucleotideType, fill = NucleotideType)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.4) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  theme_classic() +
  scale_color_manual(values = c("#56B4E9", "#D55E00", "#000000"))+
  scale_fill_manual(values = c("#56B4E9", "#D55E00", "#000000"))+
  #scale_color_colorblind() +
  #scale_fill_colorblind() +
  scale_x_discrete(labels = c("2021EE30" = "20",
                              "2021EE31" = "30",
                              "2021EE32" = "50",
                              "2021EE33" = "70",
                              "2021EE34" = "100",
                              "2021EE35" = "140",
                              "2021EE36" = "200"))+
  labs(x = "\nGeneration", y = "Nucleotide Frequency",
       title = "Frequency of SNP C > A in UMAG_05545 in the Line C"
       ) +
  #scale_fill_manual(values = c("#3C8200", "gray", "#7FAFD2")) +
  #scale_color_manual(values = c("#3C8200", "gray", "#7FAFD2")) +
  theme(axis.text.x = element_text(vjust = 0.5, color = "black", size = 10),
        axis.text.y = element_text(size = 10, color = "black"),
        plot.title = element_text(hjust = 0.30, size = 12, color = "black"),
        #plot.subtitle = element_text(hjust = 0.5),
        legend.justification = "top",
        axis.title = element_text(size = 12, color = "black"),
        #legend.position = "top" ,
        legend.text = element_text(size = 8)
        ) +
  guides(fill = guide_legend(title = "Nucleotide"),
         color = guide_legend(title = "Nucleotide",
                              title.theme = element_text(size = 9)))
  
# Concatenate both figures into a single panel
Figure.Variants.Pools <- plot_grid(Figure.SNPs.Pools, Figure.Chr9.Pools, ncol = 2, labels = c("A)", "B)") )
# Figure4 <- plot_grid(Figure.4A, Figure.4B, ncol = 2, labels = c("A)", "B)") )

plot.Figure.Variants.Pools <- "Figure_Variants.Pools.png"

ggsave(filename = plot.Figure.Variants.Pools, plot = Figure.Variants.Pools, width = 10.5, height = 6, units = "in", dpi = 300)