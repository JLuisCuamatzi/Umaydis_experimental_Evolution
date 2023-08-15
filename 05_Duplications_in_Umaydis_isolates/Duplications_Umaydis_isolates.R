# created by j cuamatzi
## script for supplementary figure 5

# load libraries
libraries <- c("data.table", "dplyr", "ggplot2", "tidyr", "scales")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

# Clean Env
rm(list = ls())

# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# 
umaydis.isolates <- fread("USMA_PopGen_ID.csv")
umaydis.isolates$IsolateName <- paste0(umaydis.isolates$`#ID`, "_", umaydis.isolates$Sample_name )

# for to read the files
for (isolate in umaydis.isolates$IsolateName){
  # create file name
  file2read <- paste0("normalizedCoverage/", isolate, ".NormCov.txt")
  # Name for the df
  df.name <- paste0("df.", isolate)
  # Read file
  df.tmp <- fread(file2read)
  # Add column with sample
  df.tmp$Sample <- isolate
  
  # assign as df
  assign(df.name, df.tmp)
  
  rm(df.name, df.tmp)
}

rm(file2read, isolate)

umaydis.pop.df <- bind_rows(mget(ls(pattern = "df.")))

rm(list = ls(pattern = "df."))




# remove the mt (24) and non nuclear contigs (25..28)
umaydis.pop.df <- umaydis.pop.df[!chr %in% c("USMA_521_v2_24", "USMA_521_v2_25", "USMA_521_v2_26", "USMA_521_v2_27", "USMA_521_v2_28")]

# Create an ID for each sample
umaydis.pop.df$ID <- gsub("_USMA_PopGen_JYD", "", umaydis.pop.df$Sample) 
umaydis.pop.df$ID <- gsub("_USMA_PopGen_GD", "", umaydis.pop.df$ID)


umaydis.pop.df$ID <- factor(umaydis.pop.df$ID, levels = c(
  "A_Irapuato", "B_Irapuato", "C_Irapuato", "D_Irapuato", "E_Irapuato", "F_Oaxaca",
   "G_Oaxaca", "H_Oaxaca", "I_Oaxaca", "J_Oaxaca", "K_Pachuca", "L_Pachuca", "M_Pachuca", "N_Pachuca", "O_Pachuca",
   "P_Sinaloa", "Q_Sinaloa", "R_Sinaloa", "S_Toluca", "T_Toluca", "U_Toluca", "V_Toluca", 
  "Um_45", "Um_48", "Um_52", "Um_55", "Um_57", "Um_58",
  
  "Um_100", "Um_101", "Um_107", "Um_110", "Um_111", "Um_112", "Um_113", "Um_114", "Um_120", "Um_121", "Um_128", "Um_130", "Um_131", "Um_134", "Um_135", "Um_142", "Um_143", "Um_149", "Um_150", "Um_157", "Um_158", "Um_159"
   )
  )

# Reverse the levels
umaydis.pop.df$ID <- fct_rev(umaydis.pop.df$ID)

# modify the chr (remove 'USMA_521_v2_')
umaydis.pop.df$chr <- as.numeric(gsub("USMA_521_v2_", "", umaydis.pop.df$chr))


# Plotting as heatmap
plot.heatmap.windows <- ggplot(umaydis.pop.df, aes(y = ID, x = window.end, fill = normalized.coverage))+
  geom_tile()+
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "both")+
  theme_classic()+
  geom_hline(yintercept = seq(1.5, 50.5, 1), color = "white", linewidth = 2)+
  labs(y = "Sample\n", x = "\nChromosome") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(0, "cm"),
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(),
    
    strip.background = element_blank(),
    strip.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    #axis.title.x = element_blank(),
    axis.line = element_blank())+
  scale_fill_gradient(low="white", high="blue", limits = c(0.7,2.1)) + # GrB  
  #scale_fill_gradient(low="white", high="red", limits = c(0,2)) + # GrR
  #scale_fill_gradient(low="black", high="red", limits = c(0,2)) + # Gr
  #scale_fill_viridis(discrete=FALSE) + # V
  # without scale fill = default version
  guides(fill = guide_colorbar(title = "Normalized \nCoverage", 
                               barwidth = 0.5, ticks = T));plot.heatmap.windows


#### export plot
plot.usma.chr.dup <- paste("Umaydis.chr.duplications.50isolates.png", sep = "/")

ggsave(filename = plot.usma.chr.dup, 
       plot = plot.heatmap.windows, 
       width = 15, 
       height = 12, 
       units = "in", 
       dpi = 300)

rm(list = ls())



