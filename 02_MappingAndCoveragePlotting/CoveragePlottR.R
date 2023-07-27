# script to plot the coverage
#### script for R version = r/3.6.1
# Libraries
library_names <- c("data.table", "zoo", "ggplot2", "tidyverse", "dplyr", "optparse", "assertthat")

# Load each library. Install it requires it
for (lib in library_names) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

## Arguments
option_list <- list(
  #
  make_option(c("--normalizedCov_file"), type = "character", default = NULL, help = "File with normalized coverage"),
  make_option(c("--window_size"), type = "integer", default = 1000, help = "Numeric value to define the window size" ),
  make_option(c("--sample"), type = "character", default = NULL, help = "Name of the sample"),
  make_option(c("--chr"), type = "character", default = NULL, help = "Chromosome prefix given the reference genome"),
  # to print help if arguments are empty
  make_option(c("--verbose"), action = "store_true", default = FALSE, help = "Print extra information")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the help option is indicated or if arguments are empty
if (is.null(opt$normalizedCov_file) || is.null(opt$sample) ) {
  print_help(opt_parser)
  quit(status = 0)
}

# Set number of threads for data.table library
if( getDTthreads() < 10 || getDTthreads() > 10){
  setDTthreads(threads = 10)
}else{print(paste("WARNING: data.table package is running with ", getDTthreads(), " threads.", sep=''))}


# Read file with normalized coverage
depth2plot <- fread(opt$normalizedCov_file)
depth2plot$chr <- as.numeric(gsub("USMA_521_v2_", "", depth2plot$chr))

# Create vector for colors
colorsForPlot <- depth2plot %>% summarise(chr = unique(chr)) %>% arrange(chr) ## identify the unique chromosomes
colorsForPlot <- ifelse(colorsForPlot$chr %% 2 == 1, "black", "gray") ## create colors for plot 

# Extract median coverage (global) for this sample
sample.median <- unique(depth2plot$global.coverage.median)

# Plot: Raw Coverage
plot.1 <- ggplot(depth2plot, aes(x = window.end, y = window.median.coverage, colour = as.factor(chr))) +
   geom_point() +
   facet_grid(.~ chr, space = 'free_x', scales = "free_x", switch = "both" ) +
   theme_bw() +
   scale_y_continuous(limits = c(0, depth2plot$global.coverage.median*5)) +
   geom_hline(yintercept = unique(depth2plot$global.coverage.median), linetype = "dashed", color = "white", size = 2) + # add a line in the global media coverage
   labs(title = paste("Sample: ", opt$sample, " (", (opt$window_size/1000), " non-overlapping kb)", sep = ""),
        subtitle = paste0("Global Median Coverage = ", sample.median),
        x = "\nChromosome",
        y = "Median Coverage Depth\n") +
   scale_colour_manual(values = colorsForPlot) +
   theme(legend.position = "none",
         panel.spacing.x = grid::unit(0, "cm"),
         panel.border = element_rect(colour = "grey", size = 0.1), panel.ontop = FALSE,
         axis.title.x = element_text(face = "bold", color = "black", size = 22),
         axis.title.y = element_text(face = "bold", color = "black", size = 22),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text.y = element_text(size = 16, color = "black"),
         plot.title = element_text(face = "bold", color = "darkred", size = 24, hjust = 0.5),
         plot.subtitle = element_text(hjust = 0.5, size = 18, color = "black"),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         strip.placement = "bottom",
         strip.text = element_text(size = 18, face = "bold")) # make the plot

# Plot: Normalized Coverage
plot.2 <- ggplot(depth2plot, 
                     aes(x = window.end, y = normalized.coverage, color = as.factor(chr))) +
  geom_point(alpha = 0.4) +
  facet_grid(.~ chr, scales = "free", space = "free", switch = "both") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "white") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "darkgray") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "darkgray") +
  geom_hline(yintercept = 4, linetype = "dashed", color = "darkgray") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgray") +
  theme_bw() +
  labs(title = paste("Sample: ", opt$sample, " (", (opt$window_size/1000), " non-overlapping kb)", sep = ""),
       subtitle = paste0("Global Median Coverage = ", sample.median),
        x = "\nChromosome",
        y = "Normalized Coverage\n") +
  scale_y_continuous(limits = c(0,5), breaks = seq(1,5,1))+
  scale_colour_manual(values = colorsForPlot) +
  theme(legend.position = "none",
         panel.spacing.x = grid::unit(0, "cm"),
         panel.border = element_rect(colour = "grey", size = 0.1), panel.ontop = FALSE,
         axis.title.x = element_text(face = "bold", color = "black", size = 22),
         axis.title.y = element_text(face = "bold", color = "black", size = 22),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text.y = element_text(size = 16, color = "black"),
         plot.title = element_text(face = "bold", color = "darkred", size = 24, hjust = 0.5),
         plot.subtitle = element_text(hjust = 0.5, size = 18, color = "black"),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         strip.placement = "bottom",
         strip.text = element_text(size = 18, face = "bold"))
  
  
# Export plots as png
# Change directory
setwd(dirname(opt$normalizedCov_file))

# Directory to save the plots
dir2save <- "CoveragePlots/"

# Create output directory for plots
if (dir.exists(dir2save)){
  print("Directory for Plots Already Exists!")
} else {
  dir.create(dir2save)
}

# Raw
plot2save.raw <- paste0(dir2save, opt$sample, "_RawCoverage.", (opt$window_size/1000),  "kb.png")
ggsave(plot2save.raw , plot = plot.1, width = 20, height = 12, units = "in", dpi = 300) # save the plot

# Normalized
plot2save.normalized <- paste0(dir2save, opt$sample, "_NormalizedCoverage.", (opt$window_size/1000),  "kb.png")
ggsave(plot2save.normalized, plot = plot.2, width = 20, height = 12, units = "in", dpi = 300) # save the plot

rm(list = ls())

