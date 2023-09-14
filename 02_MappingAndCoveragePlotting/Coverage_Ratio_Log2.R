# Script for Figure 3
# author: jcuamatzi

# load libraries
libraries <- c("data.table", "dplyr", "tidyr", "ggplot2", "scales", "cowplot", "stringr")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}


rm(list = ls()) # clean env.

# set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


### Functions
# Function to read the files with normalized coverage data:

read.NormCov.file <- function(sample){
  file2read <- paste0(sample, "_NormalizedCoverage.txt")
  
  df.tmp <- fread(file2read)
  df.tmp$Sample <- paste0(sample)
  
  return(df.tmp)
  
}

# Function to obtained the adjusted normalized coverage:

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

# Function to compute log2
calculate_log2ratio <- function(data.set, ref_col, end_col) {
  
  data.set <- data.set %>%
    # Select the next columns:
    select(Sample, chr, window.start, window.end, Adj.Norm.Cov) %>% 
    # Transform the table with pivot_wider
    pivot_wider(names_from = Sample, values_from = Adj.Norm.Cov) %>% 
    # Divide each value in the indicated columns / reference column (SG200)
    mutate(across(c(!!sym(ref_col): !!sym(end_col)),
                  .fns = ~./ !!sym(ref_col),
                  .names = "{.col}_Ratio" ))
  
  ##
  ##
  ref_col_Ratio <- paste0(ref_col, "_Ratio")
  end_col_Ratio <- paste0(end_col, "_Ratio")
  
  data.set <- data.set %>% 
    select(chr, window.start, window.end, !!sym(ref_col_Ratio):!!sym(end_col_Ratio)) %>% 
    pivot_longer(cols = !!sym(ref_col_Ratio):!!sym(end_col_Ratio), 
                 names_to = "sample", values_to = "ratio") %>% 
    mutate(log2ratio = log2(ratio))
  
  data.set$sample <- gsub("_Ratio", "", data.set$sample)
  
  return(data.set)
  
}

## Function to plot the log2 ratio as heatmap
plot.heatmap.log2 <- function(data.set, target.chr){
  # Subset the data set to a specific chr
  df.AllSamples2plot.1 <- data.set %>% 
    filter(chr == target.chr)
  
  # Create a label for the x-axis 
  chr.label<-unique( paste0("\n",gsub("USMA_521_v2_", "Chromosome ", df.AllSamples2plot.1$chr), " (kb)"))
  
  # Plot with ggplot
  my.plot <- ggplot(df.AllSamples2plot.1) +
    geom_tile(aes(x = (window.end/1000), 
                  y = y.order, 
                  fill = log2ratio)) +
    #geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), color = "white", linewidth = 2)+
    facet_grid(~ Line, space = "free", scales = "free") +
    geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), color = "white", linewidth = 2)+
    facet_grid(~ Line, space = "free", scales = "free") +
    scale_x_continuous(limits = c(0, max(df.AllSamples2plot.1$window.end)/1000),
                       breaks = c(seq(0, max(df.AllSamples2plot.1$window.end)/1000, 150), 
                                  (max(df.AllSamples2plot.1$window.end)/1000))) +
    scale_fill_gradient(low="white", high="blue", limits = c(-1,3)) + # GrB  
    #scale_fill_gradient(low="white", high="blue", limits = c(0,3.5)) + # GrB  original
    theme_classic()+
    labs(x = chr.label, y = "Colonies\n") +
    theme(
      # titles
      plot.title = element_blank(),
      axis.title.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 12, color = "black"),
      # axis text
      axis.text.x = element_text(size = 9, color = "black", angle = 90),
      axis.text.y = element_text(size = 9, color = "black"),
      # axis ticks
      axis.ticks.y = element_blank(),
      
      #axis.text = element_blank(),
      #axis.text.y = element_blank(),
      axis.line = element_blank(),
      # strip
      strip.background = element_blank(),
      strip.text = element_text(size = 12, color = "black"),
      legend.position = "right") +
    guides(fill = guide_colorbar(title = expression("log"["2"]*" (Ratio)"),
                                 barwidth = 0.5, ticks = T))
  
  
  return(my.plot)
  
}




# reading files with sample information
df.samples <- fread("../USMA_EE_Colonies_SampleSheet.csv")

# Read whole info of samples
df.experiment.info <- fread("../USMA_EE_Information.csv")
# Split Sample col in two
split_data <- str_split(df.experiment.info$Sample, "_", n = 2, simplify = TRUE)
df.experiment.info$Sample <- split_data[, 1]
df.experiment.info$ColName <- split_data[, 2]
rm(split_data)

# Read files with normalized coverage

setwd("normalizedCoverageTables/")


for (sample in df.samples$SampleID){
  
  # Create object with the name of the df
  df.name <- paste0("df.", sample)
  
  # Check if the file exists before reading it
  file2read <- paste0(sample, "_NormalizedCoverage.txt")
  
  if (file.exists(file2read)) {
    # Read file into a tmp df
    df.tmp <- read.NormCov.file(sample = sample)
    
    # Assign df.name and df.tmp
    assign(df.name, df.tmp)
    
    rm(df.name, df.tmp)
    
  } else {
    
    # File does not exist, print a message and continue with the next sample
    cat("File not found:", file2read, "\n")
  }
}

# List objects based on the pattern: df.2021EE
df.list <- mget(ls(pattern = "df.2021EE"))

# Merge all the data frames in the list using rbind
df.AllSamples <- do.call(rbind, df.list)

# Remove tmp files
rm(list = ls(pattern = "df.2021EE"))

# Adjusted the normalized coverage
df.AllSamples <- adjust.coverage(dataset.coverage = df.AllSamples, threshold.cov = 0.25 )

# Compute log2 ratio 
df.AllSamples <- calculate_log2ratio(data.set = df.AllSamples, ref_col = "2021EE01", end_col = "2021EE27")

# Add sample info to df.cov.chr9
df.AllSamples <- df.AllSamples %>% 
  left_join(select(df.experiment.info, Sample, ColName, Line, Group, H2O2Shock),
            by = c("sample" = "Sample")) 

# Remove non-nuclear chr
df.AllSamples <- df.AllSamples %>% filter(!chr %in% c("USMA_521_v2_24","USMA_521_v2_25",
                                                      "USMA_521_v2_26","USMA_521_v2_27",
                                                      "USMA_521_v2_28"))

#
#


## Remove samples: 2021EE02, 2021EE03, 2021EE04 (colonies from gen 100) and SG200
# The plots will only show the data of colonies at 200 gen.
df.AllSamples2plot <- df.AllSamples %>% filter(!sample %in% c("2021EE02", "2021EE03", "2021EE04", "2021EE01", "2021EE23", "2021EE24"))

## Indicate the y.order for the plot
df.AllSamples2plot <- df.AllSamples2plot %>% 
  mutate(y.order = case_when(
    endsWith(df.AllSamples2plot$ColName, ".C1") ~ "1",
    endsWith(df.AllSamples2plot$ColName, ".C2") ~ "2",
    endsWith(df.AllSamples2plot$ColName, ".C3") ~ "3",
    endsWith(df.AllSamples2plot$ColName, ".C4") ~ "4",
    endsWith(df.AllSamples2plot$ColName, ".C5") ~ "5"
  ))

df.AllSamples2plot$y.order <- factor(df.AllSamples2plot$y.order, levels = c("5", "4", "3", "2", "1"))
df.AllSamples2plot$Line <- gsub("Line", "Line ", df.AllSamples2plot$Line)

usma.chr <- unique(df.AllSamples$chr)


setwd("../")
dir.create("Figures_RatioLog2")
setwd("Figures_RatioLog2")

for (chr in usma.chr){
  
  plot.file <- paste0("Plot_", chr, ".Log2Ratio.png")
  
  print(paste0("Plotting: ", plot.file))
  
  my.plot <- plot.heatmap.log2(data.set = df.AllSamples2plot , target.chr = chr)
  
  ggsave(filename = plot.file, plot = my.plot, width = 8, height = 2.5, units = "in", dpi = 300)
  
}


rm(list = ls())


















