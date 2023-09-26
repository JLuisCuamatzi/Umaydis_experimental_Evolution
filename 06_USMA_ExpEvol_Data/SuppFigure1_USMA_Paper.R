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


## H2O2
# Plot A
# extract interest line
plot.target <- data.H2O2 %>% filter(Line == "Line A") %>% mutate(strip2plot = "Line A")

# data without target line
plot.no.target <- data.H2O2 %>% filter(Line != "Line A")

plot.A <- ggplot( ) +
  # plotting background lines
  geom_line(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), linewidth = 0.1) +
  geom_point(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), shape = 18, color = "gray65", size = 3) +
  scale_color_manual(values = c("gray65", "gray65")) +
  # plotting target
  geom_line(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), color = color.la, linewidth = 1, alpha = 0.8) +
  geom_point(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), shape = 18, color = color.la, size = 3) +
  # facet
  facet_grid(rows = vars(strip2plot)) +
  theme_classic() +
  geom_hline(yintercept = 1e4,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e5,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e6,linetype='dashed', col = 'darkgray') +
  geom_hline(yintercept = 1e7,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e8,linetype='dotted', col = 'gray') +
  scale_y_log10(labels = scientific, limits = c(1e4, 1e8), breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)) +
  labs(y = "Cellular concentration (cells/mL)\n", x = "\nTime (h)") +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, color = "black"),
        #strip.background = element_blank(),
        axis.title.x = element_text(face = "bold", size = 16, color = "white"),
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(face = "bold", size = 18, color = "white"),
        axis.title.y = element_blank(),
        #axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5),
        axis.text.x = element_text(color = "white", size = 16, angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 16),
        
        axis.line.x = element_blank())


# Plot B
# extract interest line
plot.target <- data.H2O2 %>% filter(Line == "Line B") %>% mutate(strip2plot = "Line B")

# data without target line
plot.no.target <- data.H2O2 %>% filter(Line != "Line B")

plot.B <- ggplot( ) +
  # plotting background lines
  geom_line(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), linewidth = 0.1) +
  geom_point(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), shape = 18, color = "gray65", size = 3) +
  scale_color_manual(values = c("gray65", "gray65")) +
  # plotting target
  geom_line(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), color = color.lb, linewidth = 1, alpha = 0.8) +
  geom_point(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), shape = 18, color = color.lb, size = 3) +
  # facet
  facet_grid(rows = vars(strip2plot)) +
  theme_classic() +
  geom_hline(yintercept = 1e4,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e5,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e6,linetype='dashed', col = 'darkgray') +
  geom_hline(yintercept = 1e7,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e8,linetype='dotted', col = 'gray') +
  scale_y_log10(labels = scientific, limits = c(1e4, 1e8), breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)) +
  labs(y = "Cellular concentration (cells/mL)\n", x = "\nTime (h)") +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, color = "black"),
        #strip.background = element_blank(),
        axis.title.x = element_text(face = "bold", size = 16, color = "white"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 16, color = "black"),
        #axis.title.y = element_blank(),
        #axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5),
        axis.text.x = element_text(color = "white", size = 16, angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 16),
        
        axis.line.x = element_blank())
  

# Plot C
# extract interest line
plot.target <- data.H2O2 %>% filter(Line == "Line C") %>% mutate(strip2plot = "Line C")

# data without target line
plot.no.target <- data.H2O2 %>% filter(Line != "Line C")

plot.C <- ggplot( ) +
  # plotting background lines
  geom_line(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), linewidth = 0.1) +
  geom_point(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), shape = 18, color = "gray65", size = 3) +
  scale_color_manual(values = c("gray65", "gray65")) +
  # plotting target
  geom_line(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), color = color.lc, linewidth = 1, alpha = 0.8) +
  geom_point(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), shape = 18, color = color.lc, size = 3) +
  # facet
  facet_grid(rows = vars(strip2plot)) +
  theme_classic() +
  geom_hline(yintercept = 1e4,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e5,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e6,linetype='dashed', col = 'darkgray') +
  geom_hline(yintercept = 1e7,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e8,linetype='dotted', col = 'gray') +
  scale_y_log10(labels = scientific, limits = c(1e4, 1e8), breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)) +
  scale_x_continuous(limits = c(0,1000), breaks = seq(0, 1000, 48)) +
  labs(y = "Cellular concentration (cells/mL)\n", x = "\nTime (h)") +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, color = "black"),
        #strip.background = element_blank(),
        axis.title.x = element_text(face = "bold", size = 16, color = "black"),
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(face = "bold", size = 18, color = "white"),
        axis.title.y = element_blank(),
        #axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5),
        axis.text.x = element_text(color = "black", size = 16, angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 16),
        
        #axis.line.x = element_blank()
        )

plot.H2O2 <- plot_grid(plot.A, plot.B, plot.C, nrow = 3, align = "v")



# Controles
# Plot D
# extract interest line
plot.target <- data.control %>% filter(Line == "Line D") %>% mutate(strip2plot = "Line D")

# data without target line
plot.no.target <- data.control %>% filter(Line != "Line D")

plot.D <- ggplot( ) +
  # plotting background lines
  geom_line(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), linewidth = 0.1) +
  geom_point(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), shape = 18, color = "gray65", size = 3) +
  scale_color_manual(values = c("gray65", "gray65")) +
  # plotting target
  geom_line(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), color = color.ld, alpha = 0.8, linewidth = 1) +
  geom_point(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), shape = 18, color = color.ld, size = 3) +
  # facet
  facet_grid(rows = vars(strip2plot)) +
  theme_classic() +
  geom_hline(yintercept = 1e4,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e5,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e6,linetype='dashed', col = 'darkgray') +
  geom_hline(yintercept = 1e7,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e8,linetype='dotted', col = 'gray') +
  scale_y_log10(labels = scientific, limits = c(1e4, 1e8), breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)) +
  labs(y = "Cellular concentration (cells/mL)\n", x = "\nTime (h)") +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, color = "black"),
        #strip.background = element_blank(),
        axis.title.x = element_text(face = "bold", size = 16, color = "white"),
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(face = "bold", size = 18, color = "white"),
        axis.title.y = element_blank(),
        #axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5),
        axis.text.x = element_text(color = "white", size = 16, angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 16),
        
        axis.line.x = element_blank())


# Plot E
# extract interest line
plot.target <- data.control %>% filter(Line == "Line E") %>% mutate(strip2plot = "Line E")

# data without target line
plot.no.target <- data.control %>% filter(Line != "Line E")

plot.E <- ggplot( ) +
  # plotting background lines
  geom_line(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), linewidth = 0.1) +
  geom_point(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), shape = 18, color = "gray65", size = 3) +
  scale_color_manual(values = c("gray65", "gray65")) +
  # plotting target
  geom_line(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), color = color.le, alpha = 0.8, linewidth = 1) +
  geom_point(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), shape = 18, color = color.le, size = 3) +
  # facet
  facet_grid(rows = vars(strip2plot)) +
  theme_classic() +
  geom_hline(yintercept = 1e4,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e5,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e6,linetype='dashed', col = 'darkgray') +
  geom_hline(yintercept = 1e7,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e8,linetype='dotted', col = 'gray') +
  scale_y_log10(labels = scientific, limits = c(1e4, 1e8), breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)) +
  labs(y = "Cellular concentration (cells/mL)\n", x = "\nTime (h)") +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, color = "black"),
        #strip.background = element_blank(),
        axis.title.x = element_text(face = "bold", size = 16, color = "white"),
        #axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 16, color = "black"),
        #axis.title.y = element_blank(),
        #axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5),
        axis.text.x = element_text(color = "white", size = 16, angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 16),
        
        axis.line.x = element_blank())


# Plot F
# extract interest line
plot.target <- data.control %>% filter(Line == "Line F") %>% mutate(strip2plot = "Line F")

# data without target line
plot.no.target <- data.control %>% filter(Line != "Line F")

plot.F <- ggplot( ) +
  # plotting background lines
  geom_line(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), linewidth = 0.1) +
  geom_point(data = plot.no.target, aes(x = cumulative_time, y = Cellular_concentration, color = Line), shape = 18, color = "gray65", size = 3) +
  scale_color_manual(values = c("gray65", "gray65")) +
  # plotting target
  geom_line(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), color = color.lf, alpha = 0.8, linewidth = 1) +
  geom_point(data = plot.target, aes(x = cumulative_time, y = Cellular_concentration), shape = 18, color = color.lf, size = 3) +
  # facet
  facet_grid(rows = vars(strip2plot)) +
  theme_classic() +
  geom_hline(yintercept = 1e4,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e5,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e6,linetype='dashed', col = 'darkgray') +
  geom_hline(yintercept = 1e7,linetype='dotted', col = 'gray') +
  geom_hline(yintercept = 1e8,linetype='dotted', col = 'gray') +
  scale_y_log10(labels = scientific, limits = c(1e4, 1e8), breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)) +
  scale_x_continuous(limits = c(0,1000), breaks = seq(0, 1000, 48)) +
  labs(y = "Cellular concentration (cells/mL)\n", x = "\nTime (h)") +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, color = "black"),
        #strip.background = element_blank(),
        axis.title.x = element_text(face = "bold", size = 16, color = "black"),
        #axis.title.x = element_blank(),
        #axis.title.y = element_text(face = "bold", size = 18, color = "white"),
        axis.title.y = element_blank(),
        #axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5),
        axis.text.x = element_text(color = "black", size = 16, angle = 90, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 16),
        
        #axis.line.x = element_blank()
  )


plot.control <- plot_grid(plot.D, plot.E , plot.F, nrow = 3, align = "v")




# 
# 
# #
# plot.H2O2 <- ggplot(data.H2O2, aes(x = cumulative_time, y = Cellular_concentration)) +
#   geom_line(aes(colour = Line)) + 
#   geom_hline(yintercept = 1e4,linetype='dotted', col = 'gray') +
#   geom_hline(yintercept = 1e5,linetype='dotted', col = 'gray') +
#   geom_hline(yintercept = 1e6,linetype='dashed', col = 'darkgray') +
#   geom_hline(yintercept = 1e7,linetype='dotted', col = 'gray') +
#   geom_hline(yintercept = 1e8,linetype='dotted', col = 'gray') +
#   scale_y_log10(labels = scientific, 
#                 limits = c(1e4, 1e8),
#                 breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)) +
#   scale_x_continuous(limits = c(0,1000), breaks = seq(0, 1000, 48)) +
#   geom_point(aes(x = cumulative_time, 
#                  y = Cellular_concentration,
#                  shape = Quantification,
#                  colour = Line))+
#   facet_grid(Line~.) +
#   theme_classic() + 
#   labs(y = "Cellular concentration (cells/mL)\n", x = "\nTime (h)") + 
#   theme(legend.position = "none",
#         strip.text = element_text(size = 12, color = "black"),
#         #strip.background = element_blank(),
#         axis.title.x = element_text(face = "bold", size = 12),
#         axis.title.y = element_text(face = "bold",
#                                     size = 12),
#         axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5),
#         axis.text.y = element_text(color = "black", size = 10),
#         
#         panel.spacing = unit(0.7, "cm")
#         )+
#   scale_shape_manual(values = c(16,2,25)) +
#   scale_colour_manual(values = c(color.la, color.lb, color.lc))
# 
# plot.H2O2
# ## control
# plot.control <- ggplot(data.control, aes(x = cumulative_time, y = Cellular_concentration)) +
#   geom_line(colour = "#0070C0") + 
#   geom_hline(yintercept = 1e4,linetype='dotted', col = 'gray') +
#   geom_hline(yintercept = 1e5,linetype='dotted', col = 'blue') +
#   geom_hline(yintercept = 1e6,linetype='dotted', col = 'red') +
#   geom_hline(yintercept = 1e7,linetype='dotted', col = 'orange') +
#   geom_hline(yintercept = 1e8,linetype='dotted', col = 'gray') +
#   scale_y_log10(labels = scientific, 
#                 limits = c(1e4, 1e8),
#                 breaks = c(1e4, 1e5, 1e6, 1e7, 1e8)) +
#   scale_x_continuous(limits = c(0,1000), 
#                      breaks = seq(0,1000,48)) +
#   geom_point(aes(x = cumulative_time, 
#                  y = Cellular_concentration,
#                  shape = Quantification), colour = "#0070C0")+
#   facet_grid(Line~.) +
#   theme_classic() + 
#   labs(y = "Cellular concentration (cells/mL)\n", x = "\nTime (h)") + 
#   theme(legend.position = "none",
#         strip.text = element_text(size = 12, color = "black"),
#         #strip.background = element_blank(),
#         axis.title.x = element_text(face = "bold", size = 12),
#         axis.title.y = element_text(face = "bold",
#                                     size = 12),
#         axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5),
#         axis.text.y = element_text(color = "black", size = 10),
#         
#         panel.spacing = unit(0.7, "cm")
#   )+
#   scale_shape_manual(values = c(16,2,25)) +
#   scale_colour_manual(values = c(color.ld, color.le, color.lf))




plot.SuppFigure1 <- plot_grid(plot.H2O2,  plot.control, 
                              nrow = 2, scale = 0.9, align = "h", 
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
       width = 13, height = 16, units = "in", dpi = 300)

rm(list = ls())







