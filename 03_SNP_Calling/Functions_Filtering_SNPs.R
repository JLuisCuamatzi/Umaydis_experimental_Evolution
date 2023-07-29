
#load libraries
library(vcfR) # works
library(tidyr) # works
library(tidyverse) # works
library(data.table)# works
library(dplyr)  # works
library(UpSetR) # works
library(ComplexHeatmap) # works
library(ggplot2) # works
library(Cairo) # works
library(cowplot) # works
library(ggthemes) # works
library(circlize)


### my R functions

# getSNP <- function (vcf.file){
#   vcf <- read.vcfR(vcf.file) # read vcf (vcfR::read.vcfR)
#   
#   vcf <- as.data.frame(bind_cols(getFIX(vcf), extract.gt(vcf))) # extract the data from vcf (chr, pos, ref, alt, sampleGT)
#   
#   vcf <- vcf[, -c(3,7)] # remove empty cols
#   
#   vcf <- vcf %>% mutate(ToRemoveRef = if_else( (nchar(vcf$REF) != 1), 0, 1 ),
#                         ToRemoveAlt = if_else( (nchar(vcf$ALT) != 1), 0, 1 )) 
#   
#   ## vcf have info of SNP and InDel. Keep only SNP based on the length of string in REF and ALT.
#   
#   vcf <- vcf %>% mutate(Rows2Remove = if_else(((ToRemoveRef + ToRemoveAlt) == 2), 1, 0))
#   
#   vcf <- vcf %>% filter(Rows2Remove == 1) %>% select(1:6) %>% pivot_longer(cols = 6, names_to = "Sample", values_to = "GT")    
#   
#   vcf$SNP_ID <- paste(vcf$CHROM, 
#                       vcf$POS,
#                       vcf$REF,
#                       vcf$ALT,
#                       vcf$GT, sep = "-")
#   
#   return(vcf)
# }


## SNP Matrix - Raw Numbers Symmetric
intersect.snp.raw <- function(list_snp){
  mtrx.comb <- combn(names(list_snp), 2)              # Create all possible combinations
  
  #mtx <- mat.or.vec( nr = length(unique(list_snp)), 
  #                   nc = length(unique(list_snp)))   # Create an empty matrix
  # change unique to names
  mtx <- mat.or.vec( nr = length(names(list_snp)), 
                     nc = length(names(list_snp)))   # Create an empty matrix
  colnames(mtx) <- names(list_snp)                    # col names to the matrix 
  rownames(mtx) <- names(list_snp)                    # row names to the matrix
  
  #mtx.perc.col <- mtx # perc (/col)
  #mtx.perc.row <- mtx # perc (/row)
  
  for (i in 1:ncol(mtrx.comb)){
    x.lab <- as.character(mtrx.comb[1,i])
    y.lab <- as.character(mtrx.comb[2,i])
    x.index <- which(names(list_snp) == x.lab)
    y.index <- which(names(list_snp) == y.lab)
    #f <- length(intersect(list_snp[[x.index]], list_snp[[y.index]]))
    # symmetric
    mtx[x.index, x.index] <- length(intersect(list_snp[[x.index]], list_snp[[x.index]]))
    mtx[y.index, y.index] <- length(intersect(list_snp[[y.index]], list_snp[[y.index]]))
    mtx[x.index, y.index] <- length(intersect(list_snp[[x.index]], list_snp[[y.index]]))
    mtx[y.index, x.index] <- length(intersect(list_snp[[y.index]], list_snp[[x.index]]))
    #
    
    # non symmetric
    # mtx[x.index, x.index] <- length(intersect(list_snp[[x.index]], list_snp[[x.index]]))
    # mtx[y.index, y.index] <- length(intersect(list_snp[[y.index]], list_snp[[y.index]]))
    # mtx[x.index, y.index] <- length(intersect(list_snp[[y.index]], list_snp[[x.index]]))
    # mtx[y.index, x.index] <- length(intersect(list_snp[[y.index]], list_snp[[x.index]]))
  }
  return(mtx)
}

## SNP Matrix - Raw Numbers Symmetric
intersect.snp.raw.NoSymm <- function(list_snp){
  mtrx.comb <- combn(names(list_snp), 2)              # Create all possible combinations
  
  #mtx <- mat.or.vec( nr = length(unique(list_snp)), 
  #                   nc = length(unique(list_snp)))   # Create an empty matrix
  # change unique to names
  mtx <- mat.or.vec( nr = length(names(list_snp)), 
                     nc = length(names(list_snp)))   # Create an empty matrix
  colnames(mtx) <- names(list_snp)                    # col names to the matrix 
  rownames(mtx) <- names(list_snp)                    # row names to the matrix
  
  #mtx.perc.col <- mtx # perc (/col)
  #mtx.perc.row <- mtx # perc (/row)
  
  for (i in 1:ncol(mtrx.comb)){
    x.lab <- as.character(mtrx.comb[1,i])
    y.lab <- as.character(mtrx.comb[2,i])
    x.index <- which(names(list_snp) == x.lab)
    y.index <- which(names(list_snp) == y.lab)
    #f <- length(intersect(list_snp[[x.index]], list_snp[[y.index]]))
    
    # non symmetric
    mtx[x.index, x.index] <- length(intersect(list_snp[[x.index]], list_snp[[x.index]]))
    mtx[y.index, y.index] <- length(intersect(list_snp[[y.index]], list_snp[[y.index]]))
    mtx[x.index, y.index] <- length(intersect(list_snp[[y.index]], list_snp[[x.index]]))
    mtx[y.index, x.index] <- length(intersect(list_snp[[y.index]], list_snp[[x.index]]))
  }
  return(mtx)
}

## SNP Matrix - No symmetric
intersect.snp.perc.NoSym <- function(list_snp){
  mtrx.comb <- combn(names(list_snp), 2)
  #mtx <- mat.or.vec( nr = length(unique(list_snp)), 
  #                   nc = length(unique(list_snp))) # raw num
  mtx <- mat.or.vec( nr = length(names(list_snp)), 
                     nc = length(names(list_snp))) # raw num
  colnames(mtx) <- names(list_snp)
  rownames(mtx) <- names(list_snp)
  for (i in 1:ncol(mtrx.comb)){
    x.lab <- as.character(mtrx.comb[1,i])
    y.lab <- as.character(mtrx.comb[2,i])
    x.index <- which(names(list_snp) == x.lab)
    y.index <- which(names(list_snp) == y.lab)
    mtx[x.index, x.index] <- length(intersect(list_snp[[x.index]], list_snp[[x.index]]))/
      length(list_snp[[x.index]])
    mtx[y.index, y.index] <- length(intersect(list_snp[[y.index]], list_snp[[y.index]]))/
      length(list_snp[[y.index]])
    mtx[x.index, y.index] <- length(intersect(list_snp[[y.index]], list_snp[[x.index]]))/
      length(list_snp[[y.index]])
    mtx[y.index, x.index] <- length(intersect(list_snp[[y.index]], list_snp[[x.index]]))/
      length(list_snp[[x.index]])
  }
  return(mtx)
}

## SNP Matrix - Symmetric
intersect.snp.perc.Sym <- function(list_snp){
  mtrx.comb <- combn(names(list_snp), 2)
  mtx <- mat.or.vec( nr = length(unique(list_snp)), 
                     nc = length(unique(list_snp))) 
  colnames(mtx) <- names(list_snp)
  rownames(mtx) <- names(list_snp)
  for (i in 1:ncol(mtrx.comb)){
    x.lab <- as.character(mtrx.comb[1,i])
    y.lab <- as.character(mtrx.comb[2,i])
    x.index <- which(names(list_snp) == x.lab)
    y.index <- which(names(list_snp) == y.lab)
    mtx[x.index, x.index] <- length(intersect(list_snp[[x.index]], 
                                              list_snp[[x.index]]))/length(list_snp[[x.index]])
    mtx[y.index, y.index] <- length(intersect(list_snp[[y.index]], 
                                              list_snp[[y.index]]))/length(list_snp[[y.index]])
    mtx[x.index, y.index] <- length(intersect(list_snp[[x.index]], 
                                              list_snp[[y.index]]))/length(list_snp[[y.index]])
    mtx[y.index, x.index] <- length(intersect(list_snp[[y.index]],
                                              list_snp[[x.index]]))/length(list_snp[[y.index]])
  }
  return(mtx)
}

## function to remove background

remove.background <- function(background.sample.df, variants.df){
  # background sample = sample to used as background
  # variants.df = a concatenated data frame with variants from all samples (SNP)
  
  setDT(variants.df) # set as data table format 
  
  # 1) extract the variants in the background sample from 'variants.df'
  df.bg <- background.sample.df
  
  # 2) Remove the background sample from `variants.df`
  #df.tmp <- variants.df[!Sample == background.sample]
  
  # 3) Filter the `variants.df`. Remove all rows that match SNP_ID in the background
  df.No.BG <- variants.df[! (variants.df$SNP_ID %in% df.bg$SNP_ID), ]
  
  
  return(df.No.BG)
  
}

# write function to know the quality distribution in a df of variants
# data frame should be have a column call QUAL
# filter.step should be a string ("Line C (with SG200))
# plot.title
proportion.plot <- function(vcf.df, qual.threshold, filter.step, plot.title){
  vcf.df$QUAL <- as.numeric(vcf.df$QUAL)
  #count the number of obsv below and above threshold
  below <- sum(vcf.df$QUAL < qual.threshold)
  above <- sum(vcf.df$QUAL > qual.threshold)
  
  # Calculate proportion in percentage
  total <- nrow(vcf.df)
  prop_below <- (below / total) * 100
  prop_above <- (above / total) * 100
  
  # Create a data frame for plotting
  prop_df <- data.frame(
  Analysis = rep(filter.step, 2),
  Threshold = c(paste("Below ", qual.threshold, sep = ""), 
                paste("Above ", qual.threshold, sep = "")),
  Proportion = c(prop_below, prop_above))
    
  # plotting percentage
  plot.threshold.1 <- ggplot(prop_df, aes(x = Analysis, 
                                        y = Proportion, fill = Threshold)) +
    geom_bar(stat = "identity") +
    labs(title = plot.title,
       x = "Step",
       y = "Proportion (%)") +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          legend.position = "none",
          axis.title.y = element_text(face = "bold"))
    
    # plotting raw numbers
  raw_df <- data.frame( Analysis = rep(filter.step, 2),
                       Threshold = c(paste("Below ", qual.threshold, sep = ""), paste("Above ", qual.threshold, sep = "")),
                       Numbers = c(below, above))
  
  plot.threshold.2 <- ggplot(raw_df, aes(x = Analysis, y = Numbers, fill = Threshold)) +
    geom_bar(stat = "identity") +
    labs(title = plot.title,
       x = "Step",
       y = "Number of SNPs") +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold"))
    
    plot.both <- plot_grid(plot.threshold.1, plot.threshold.2, scale = 0.9,
                            rel_widths = c(0.7, 1))
    
  return(plot.both)
}

#### FUNCTION TO PLOT THE QUALITY DISTRIBUTION OF SNPs

plot.quality <- function (vcf.df, plot.title, plot.subtitle, qual.threshold){
  vcf.df$QUAL <- as.numeric(vcf.df$QUAL)
  
  SNPs.num <- nrow(vcf.df)
  
  plot.distr.qual <- vcf.df %>% ggplot(aes(x = as.numeric(QUAL) )) +
    geom_density(alpha=.2, fill="#FF6666") +
    geom_vline(xintercept = as.numeric(qual.threshold), linetype = "dashed", color = "red") +
    #scale_x_continuous( limits = c(0 , max(vcf.df$QUAL) )) +
    theme_classic() +
    labs(x = "Quality", y = "Relative Abudance",
         title = plot.title,
         subtitle = paste( plot.subtitle, "\nNumber of SNPs = ", SNPs.num, sep = ""  ) )+
    theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
  
  return(plot.distr.qual)
}


# matrix raw = plot 
plot.matrix.raw <- function(matrix.raw){
  # create color function
  col_fun <- colorRamp2(c(0, (max(matrix.raw)/2), (max(matrix.raw))),
                        c("yellowgreen", "seagreen4", "navy"))
  # create groups to split the matrix (col and row)
  column_split <- c(rep("A", 6), rep("B", 6), rep("C", 8),
                     rep("C.NH", 5), rep("C.H", 5),
                     rep("D", 2), rep("E", 2), rep("F", 2))
  column_split <- factor(column_split, 
                          levels = c("A", "B", "C", "C.NH", "C.H", "D", "E", "F"))
  # column_split <- c(rep("A", 6), rep("B", 6), rep("C", 8),
  #                   
  #                   rep("D", 2), rep("E", 2), rep("F", 2))
  # column_split <- factor(column_split, 
  #                        levels = c("A", "B", "C",  "D", "E", "F"))
  
  
  row_split <- c(rep("A", 6), rep("B", 6), rep("C", 8),
                 rep("C.NH", 5), rep("C.H", 5),
                 rep("D", 2), rep("E", 2), rep("F", 2))
  row_split <- factor(row_split, 
                      levels = c("A", "B", "C", "C.NH", "C.H", "D", "E", "F"))
  
  # row_split <- c(rep("A", 6), rep("B", 6), rep("C", 8),
  #                
  #                rep("D", 2), rep("E", 2), rep("F", 2))
  # row_split <- factor(row_split, 
  #                     levels = c("A", "B", "C", "D", "E", "F"))
  
  
  matrix.raw.plot <- Heatmap(matrix.raw, 
                             name = "Shared SNP\n(Raw Numbers)", 
                             cluster_rows = F, 
                             cluster_columns = F,
                             column_split = column_split,
                             row_split = row_split,
                             border = T,
                             row_gap = unit(2, "mm"),
                             column_gap = unit(2, "mm"),
                             col = col_fun,
                             cell_fun=function(j, i, x, y, w, h, col){
                               grid.text(matrix.raw[i, j], x, y)
                             })
  
  
  return(matrix.raw.plot)
}

# matrix percentage = plot
plot.matrix.per <- function(matrix.per){
  col_fun <- colorRamp2(c(0, 50, 100), c("yellowgreen", "seagreen4", "navy"))
  # create groups to split the matrix (col and row)
  
  column_split <- c(rep("A", 6), rep("B", 6), rep("C", 8),
                    rep("C.NH", 5), rep("C.H", 5),
                    rep("D", 2), rep("E", 2), rep("F", 2))
  column_split <- factor(column_split, 
                         levels = c("A", "B", "C", "C.NH", "C.H", "D", "E", "F"))
  
  # column_split <- c(rep("A", 6), rep("B", 6), rep("C", 8),
  #                   
  #                   rep("D", 2), rep("E", 2), rep("F", 2))
  # column_split <- factor(column_split, 
  #                        levels = c("A", "B", "C",  "D", "E", "F"))
  
  
  row_split <- c(rep("A", 6), rep("B", 6), rep("C", 8),
                 rep("C.NH", 5), rep("C.H", 5),
                 rep("D", 2), rep("E", 2), rep("F", 2))
  row_split <- factor(row_split, 
                      levels = c("A", "B", "C", "C.NH", "C.H", "D", "E", "F"))
  
  # row_split <- c(rep("A", 6), rep("B", 6), rep("C", 8),
  #                
  #                rep("D", 2), rep("E", 2), rep("F", 2))
  # row_split <- factor(row_split, 
  #                     levels = c("A", "B", "C", "D", "E", "F"))
  
  # prepare matrix for plotting
  matrix.per <- (matrix.per*100) # multiply per 100 
  
  matrix.per <- as.data.frame(matrix.per) %>% 
    mutate(across(where(is.numeric), ~ round(., 0))) # round the percentage to integer 
  
  matrix.per <- as.matrix(matrix.per) # turn into matrix 
  
  matrix.perc.plot <- Heatmap(matrix.per,
          name = "Shared SNP (%)", 
          cluster_rows = F, 
          cluster_columns = F,
          column_split = column_split,
          row_split = row_split,
          border = T,
          row_gap = unit(2, "mm"),
          column_gap = unit(2, "mm"),
          col = col_fun,
          cell_fun = function(j, i, x, y, w, h, col){
            grid.text(matrix.per[i, j], x, y)
          })
  
  return(matrix.perc.plot)
}


# this function uses vcfR library (install from bioconductor)
get.SNP <- function(vcf.file){
  vcf.obj <- read.vcfR(file = vcf.file) # read vcf file using vcfR::read.vcfR()
  vcf.df <- vcfR2tidy(vcf.obj) # convert to a tidy object (as list)
  vcf.df.fix <- vcf.df$fix # extract "fix" and convert to data frame
  vcf.df.fix <- vcf.df.fix %>% select(CHROM, POS, REF, ALT, QUAL, DP, INDEL) # select cols
  vcf.df.gt  <- vcf.df$gt # extract "gt" and convert to df
  vcf.df.gt  <- vcf.df.gt %>% select(POS, Indiv, gt_GT) # select cols
  # bind cols:
  vcf.df <- bind_cols(vcf.df.fix, vcf.df.gt) # concat dfs (fix and gt)
  vcf.df <- vcf.df[, -8] # remove duplicated col (POS)
  names(vcf.df)[2] <- "POS" # rename POS...2 to POS
  names(vcf.df)[8] <- "Sample" # rename Indiv to Sample
  names(vcf.df)[9] <- "GT" # rename gt_GT to GT
  
  rm(vcf.df.fix, vcf.df.gt) # remove tmp df 
  vcf.df <- vcf.df %>% filter(INDEL == FALSE) # keep only snps (remove INDELS)
  vcf.df <- vcf.df[,-7] # remove col of INDEL (T o F)
  
  
  vcf.df$QUAL <- as.numeric(vcf.df$QUAL) # set as numeric
  vcf.df$DP <- as.numeric(vcf.df$DP) # set as numeric
  vcf.df$GT <- as.numeric(vcf.df$GT) # set as numeric
  
  # finally create a unique ID for each SNP
  vcf.df$SNP_ID <- paste(vcf.df$CHROM, 
                         vcf.df$POS,
                         vcf.df$REF,
                         vcf.df$ALT,
                         vcf.df$GT, sep = "-")
  
  # return data frame seted as data table
  setDT(vcf.df)
  
  return(vcf.df)
}

## create function to plot QUAL and DP
# the input for this function is a vcf converted into a df
plot.qual.and.dp <- function(vcf.df, qual.threshold = 0, dp.threshold = 0,
                             plot.title.q, plot.title.dp) {
  # first, estimate the number of SNPs in df by counting number of rows
  SNPs.num <- nrow(vcf.df)
  min.q <- min(vcf.df$QUAL)
  max.q <- max(vcf.df$QUAL)
  median.q <- median(vcf.df$QUAL)
  n.samples <- length(unique(vcf.df$Sample))
  
  plot.q <- vcf.df %>% 
    ggplot() +
    geom_density(aes(x = QUAL), alpha = 0.2, fill = "#FF6666") +
    geom_vline(xintercept = qual.threshold, linetype = "dashed", color = "red") +
    theme_classic() +
    labs(x = "Quality", y = "Relative Abudance",
         title = plot.title.q,
         subtitle = paste(
           "Number of SNPs = ", SNPs.num, 
           "   Min Q = ", min.q,
           "   Max Q = ", max.q,
           "   Median Q = ", median.q,
           "   # Sample(s) = ", n.samples,
           sep = ""  ) )+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  min.dp = min(vcf.df$DP)
  max.dp = max(vcf.df$DP)
  median.dp = median(vcf.df$DP)
  
  plot.dp <- vcf.df %>% 
    ggplot() +
    geom_density(aes(x = DP), alpha = 0.2, fill = "blue") +
    geom_vline(xintercept = dp.threshold, linetype = "dashed", color = "blue") +
    theme_classic() +
    labs(x = "DP", y = "Relative Abudance",
         title = plot.title.dp,
         subtitle = paste(
           "Number of SNPs = ", SNPs.num, 
           "   Min DP = ", min.dp,
           "   Max DP = ", max.dp,
           "   Median DP = ", median.dp,
           "   # Sample(s) = ", n.samples,
           sep = ""  ) )+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # concat plots with cowplot::plot_grid()
  plot.q.dp <- plot_grid(plot.q, plot.dp, scale = 0.9)
  
  return(plot.q.dp)
  
}


# function to transform the input files with allele depth data
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
  
  #df_split$NucleotideType <- factor(df_split$NucleotideType, levels = c("Reference", "Alternative"))
  # return pivoted df
  return(df_split)
}

# function to plot the allele frequency

plotting.alleleFreq <- function(snp.target, df.all.samples){
  df.snp <- df.all.samples %>% filter(SNP_ID == snp.target)
  
  # df.snp <- df.snp %>% mutate(NucleotideType = if_else((Nucleotide == Alt_Detected), 
  #                                                      paste0("Alt (", Alt_Detected, ")") ,
  #                                                      NucleotideType))
  # 
  # df.snp <- df.snp %>% 
  #   mutate(NucleotideType = if_else(Nucleotide != Alt_Detected & Nucleotide != Ref, "Other",
  #                                   NucleotideType))
  
  # plotting
  plot.frequency <- df.snp %>% 
    ggplot(aes(x = Sample, y = PerFreq, 
               color = NucleotideType, fill = NucleotideType)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.4)+
    #facet_grid(~ SNP_ID_1, scales = "free") +
    # facet by line
    facet_grid(~Line, scale = "free", space = "free")+
    labs(title = paste0("Allele Frequency of ", unique(df.snp$SNP_ID)),
         subtitle = "All Colonies",
         x = "Sample", y = "Allele Frequency") +
    scale_y_continuous(labels = percent)+
    theme_classic() +
    scale_fill_manual(values = c("#3C8200", "gray", "#7FAFD2")) +
    scale_color_manual(values = c("#3C8200", "gray", "#7FAFD2")) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.title = element_text(size = 14, color = "black", face = "bold"))
  
  return(plot.frequency)
  
}



