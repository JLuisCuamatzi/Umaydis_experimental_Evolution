# created by j cuamatzi
# script for supp. figure 1

#libraries
libraries <- c("data.table", "dplyr", "tidyr", "ggplot2", "scales", "cowplot")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

rm(list = ls())

# Set as working directory, the directory in which this script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# function for scientific notation
scientific <- function(x){
  ifelse(x==0, "0", 
         parse(text=gsub("[+]", "",
                         gsub("e", "%*%10^", scientific_format()(x)))))
}

# colors 
#color.sg <- "#222222"

color.la <- "#843C0C"
color.lb <- "#7F6000"
color.lc <- "#C55A11"
color.ld <- "#002060"
color.le <- "#1F4E79"
color.lf <- "#0070C0"

# Read data set
df <- fread("F01_ExpEvol_CFU_Data/S1_ExpEvol_GrowthData.csv")

# create data frame for both groups (control and h2o2)
data.control <- df[Condition == "Control"]
data.H2O2 <- df[Condition == "H2O2"]



#
plot.H2O2 <- ggplot(data.H2O2, aes(x = cumulative_time, y = Cellular_concentration)) +
  geom_line(aes(colour = Line)) + 
  geom_hline(yintercept = 1e4,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e5,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e6,linetype='dashed', col = 'darkgray') +
  geom_hline(yintercept = 1e7,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e8,linetype='dotted', col = 'gray') +
  scale_y_log10(labels = scientific, 
                limits = c(1e4, 1e8),
                breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)) +
  scale_x_continuous(limits = c(0,1000), breaks = seq(0, 1000, 48)) +
  geom_point(aes(x = cumulative_time, 
                 y = Cellular_concentration,
                 shape = Quantification,
                 colour = Line))+
  facet_grid(Line~.) +
  theme_classic() + 
  labs(y = "Cellular concentration (cells/mL)\n", x = "\nTime (h)") + 
  theme(legend.position = "none",
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_blank(),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold",
                                    size = 14),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10))+
  scale_shape_manual(values = c(16,2,25)) +
  scale_colour_manual(values = c(color.la, color.lb, color.lc))


## control
plot.control <- ggplot(data.control, aes(x = cumulative_time, y = Cellular_concentration)) +
  geom_line(colour = "#0070C0") + 
  geom_hline(yintercept = 1e4,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e5,linetype='dotted', col = 'blue') +
  geom_hline(yintercept = 1e6,linetype='dotted', col = 'red') +
  geom_hline(yintercept = 1e7,linetype='dotted', col = 'orange') +
  geom_hline(yintercept = 1e8,linetype='dotted', col = 'gray') +
  scale_y_log10(labels = scientific, 
                limits = c(1e4, 1e8),
                breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)) +
  scale_x_continuous(limits = c(0,1000), 
                     breaks = seq(0,1000,48)) +
  geom_point(aes(x = cumulative_time, 
                 y = Cellular_concentration,
                 shape = Quantification), colour = "#0070C0")+
  facet_grid(Line~.) +
  theme_classic() + 
  labs(y = "Cellular concentration (cells/mL)\n", x = "\nTime (h)") + 
  theme(legend.position = "none",
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_blank(),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold",
                                    size = 14),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 10))+
  scale_shape_manual(values = c(16,2,25)) +
  scale_colour_manual(values = c(color.ld, color.le, color.lf))


plot.SuppFigure1 <- plot_grid(plot.H2O2,  plot.control, 
                              ncol = 2, scale = 0.9, align = "h", 
                              labels = c("A)", "B)"))


## save plot

dirSavePlots <- "Figures"

if ( dir.exists(dirSavePlots) ){
  print ("Directory already exists!!")
} else {
  dir.create(dirSavePlots)
}

plot.SuppFig1 <- paste(dirSavePlots, "SupplementaryFigure1_USMA_Paper.svg", sep = "/")

ggsave(filename = plot.SuppFig1, plot = plot.SuppFigure1,
       width = 12, height = 5, units = "in", dpi = 300)

rm(list = ls())







