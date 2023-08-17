libraries <- c("data.table", "ggplot2", "tidyr", "cowplot", "dplyr")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

rm(list = ls())

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## Create data frame with coordinates of genes located between the 140 kb to the 154 kb in chr 9
df.genes <- data.frame(
  Gene = c( "UMAG_10438", "UMAG_03439", "UMAG_03440", "HobS", "UMAG_10439", "UMAG_03442"),
  xstart = c(140.478, 143.381, 144.651, 147.099, 149.932, 151.640),
  xend = c(142.766, 144.134, 145.358, 147.497, 150.940, 153.355),
  ystart = rep(5.2, 6),
  yend = rep (5.5, 6))

df.genes <- df.genes %>% mutate(TextPositionX = (xstart + xend)/2, TextPositionY = ((ystart + yend) + 1 )/2)

df.genes$Gene <- gsub("_", "\n", df.genes$Gene)

# Plotting GC of chromosome 9
# GC was estimated using non-overlapping windows of 0.2 kb (200 bp)
df.Chr9.GC <- fread("../USMA_Genome/USMA_521_v2_9_GC_200bp.txt" )# windows of 200 bp
  #df.Chr9.GC <- fread("../USMA_521_v2_9_GC_100bp.txt") # windows of 50 bp
df.Chr9.GC2plot <- df.Chr9.GC %>% filter(End >= as.numeric(140000) & End <= as.numeric(155000))  

plot.Chr9.GC <- ggplot() +
  geom_line(data= df.Chr9.GC2plot, aes(x = (End/1000), y = GC ), linewidth = 0.5) +
  geom_hline(yintercept = 54.04, linetype = "dotted", color = "gray", linewidth = .9) + # mean gc for umaydis = 54.04%
  theme_classic()+
  labs(y = "\n\nGC (%)") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(face = "bold", color = "black", size = 8),
    axis.text.y = element_text(size = 8, color = "black")
  ) +
  geom_rect(aes(xmin = 147.099, xmax = 147.498, ymin = 20, ymax = 80), fill = "darkblue", alpha = 0.1 )+
  scale_y_continuous(limits = c(20, 80))


## Read files with normalized coverage at bp resolution
samples.2read <- c("2021EE01", "2021EE16", "2021EE18")



for (sample in samples.2read){
  file.2read <- paste0("normalizedCoverageTables/bpResolution/", sample, ".bpNormalizedCoverage.Chr9.txt.gz")
  
  df.name <- paste0("df.", sample)
  df.tmp <- fread(file.2read)
  df.tmp <- df.tmp[Position >= 140000 & Position <=155000]
  
  assign(df.name, df.tmp)
  
  rm(df.name, df.tmp, file.2read)
}


plot.coverage.bp <- ggplot()+
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 2, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 3, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_hline(yintercept = 4, linetype = "dashed", color = "gray", linewidth = 0.4)+
  geom_line(data = df.2021EE01, aes(x = (Position/1000), y = normalized.coverage), color = "purple", linewidth = 0.4) +
  geom_line(data = df.2021EE16, aes(x = (Position/1000), y = normalized.coverage), color = "orange", linewidth = 0.4) +
  geom_line(data = df.2021EE18, aes(x = (Position/1000), y = normalized.coverage), color = "darkred", linewidth = 0.4)+
    scale_y_continuous(limits = c(0, 6), breaks = seq(0, 5, 1)) +
  scale_x_continuous(breaks = c(140, 142.5, 145, 147.5, 150, 152.5, 155))+
    geom_rect(data = df.genes, aes(xmin = xstart, xmax = xend, 
                                   ymin = ystart, ymax = yend),
              fill = c("gray", "gray", "gray", "darkblue", "gray", "gray"), alpha = 0.5, color = "black") +
  theme_classic() +
  geom_text(data = df.genes, aes(x = TextPositionX, label = Gene, y = TextPositionY),
              size = 2.6) +
  labs(y = "\n\nNormalized Coverage", x = "Chromosome 9 (kb)") +
  theme_classic() +
  geom_rect(aes(xmin = 147.099, xmax = 147.498, ymin = 0, ymax = 5.2), fill = "darkblue", alpha = 0.1 )+ # HobS
  geom_rect(aes(xmin = 148.500, xmax = 149.100, ymin = 0, ymax = 5.2), fill = "darkred", alpha = 0.1 )+ #break point
  theme(
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(size = 8, color = "black"),
      
      # axis titles
      axis.title.x = element_text(face = "bold", size = 9, color = "black"),
      axis.title.y = element_text(face = "bold", size = 9, color = "black")  )

#plot.coverage.bp
plot.chr9.breakpoint <- plot_grid(plot.Chr9.GC, plot.coverage.bp, nrow = 2, align = "hv", rel_heights = c(0.4, 1), labels = c("A)", "B)"))  

ggsave(filename = "normalizedCoverageTables/CoveragePlots/Supp.Figure.NormCov.bpResolution.png", plot = plot.chr9.breakpoint, width = 16, height = 11, units = "cm", dpi = 300)  
  
rm(list = ls())
  
  
  
  
  
  
  
  
  
  
  
  
  
