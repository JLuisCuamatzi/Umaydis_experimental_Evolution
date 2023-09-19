# author: jcuamatzi

# Libraries:
libraries <- c("ggplot2", "data.table", "dplyr", "tidyr", "scales")


# Load each library. Install it requires it
for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}


rm(list = ls())

# Set as working directory, the directory in which this script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read data of CFU
df <- fread("F01_ExpEvol_CFU_Data/Umaydis_ExpEvol_GrowthData.csv")

# keep only CFU
df <- df[,-5]

# here we get only data for SG200 and treatment 20
df.Figure.1 <- df %>%  
  filter(TreatmentNum %in% c("SG200", "Treatment 20")) %>% filter(H2O2 %in% c("0 mM", "5 mM", "60 mM",
                                                                                   "w/o ROS 0 mM",
                                                                                   "w/o ROS 5 mM",
                                                                                   "w/o ROS 60 mM"))
# remove extrac charactets
df.Figure.1$H2O2 <- gsub("w/o ROS ", "", df.Figure.1$H2O2)

# keep only 0 mM, 5 mM and 60 mM
df.Figure.1 <- df.Figure.1 %>% 
  filter(H2O2 %in% c("0 mM", "5 mM", "60 mM"))

df.Figure.1 <- df.Figure.1 %>% mutate(CFU_Type =
                                        if_else(H2O2 != "0 mM", "CFU_Target", "CFU_Reference"))

df.Figure.1$H2O2 <- gsub(" ", "", df.Figure.1$H2O2)

# Pivot the table 
df.Figure.1 <- df.Figure.1 %>% 
  group_by(Condition, TreatmentNum, Line) %>% 
  summarise(CFU_0mM = CFU[H2O2 == "0mM"],
            CFU_5mM = CFU[H2O2 == "5mM"],
            CFU_60mM = CFU[H2O2 == "60mM"]) %>% 
  ungroup()

# estimate the percentage of cfu 
# we calculate 5 and 60 mM against 0 mM in each group 
df.Figure.1 <- df.Figure.1 %>% 
  mutate(across(c(`CFU_0mM`:`CFU_60mM`),
                .fns = ~./`CFU_0mM`,
                .names = "{.col}_Ratio"))

df.Figure.1 <- df.Figure.1 %>% 
  select(Condition, TreatmentNum, Line, CFU_0mM_Ratio: CFU_60mM_Ratio) %>% 
  pivot_longer(cols = CFU_0mM_Ratio: CFU_60mM_Ratio, names_to = "H2O2", values_to = "CFU_Ratio") %>% 
  setDT()

df.Figure.1$H2O2 <- gsub("CFU_", "", df.Figure.1$H2O2)
df.Figure.1$H2O2 <- gsub("_Ratio", "", df.Figure.1$H2O2)



# Keep just 5 mM and 60 mM
df.Figure.1 <- df.Figure.1[H2O2 %in% c("5mM", "60mM")]

df.Figure.1.Mean <- df.Figure.1 %>% 
  group_by(Condition, TreatmentNum, H2O2, Line) %>% 
  summarise(CFU_Ratio_Mean = mean(CFU_Ratio),
            CFU_Ratip_DS = sd(CFU_Ratio)) %>% 
  ungroup()

df.Figure.1 %>% 
  group_by(Condition, TreatmentNum, H2O2, Line) %>% 
  summarise(CFU_Ratio_Mean = mean(CFU_Ratio),
            CFU_Ratip_DS = sd(CFU_Ratio)) %>% 
  ungroup()


df2plot.1 <- df.Figure.1.Mean %>% 
  select(Condition, H2O2, Line, CFU_Ratio_Mean)

df2plot.1$Line <- factor(df2plot.1$Line, levels = c("SG200", "A", "B", "C", "D", "E", "F"))

df2plot.1 <- df2plot.1 %>% 
  mutate(Strip.Labs = 
           case_when(
             startsWith(as.character(df2plot.1$Line), "SG200") ~ "Initial\nStrain",
             startsWith(as.character(df2plot.1$Line), "A") ~ "Exposed Pools\n",
             startsWith(as.character(df2plot.1$Line), "B") ~ "Exposed Pools\n",
             startsWith(as.character(df2plot.1$Line), "C") ~ "Exposed Pools\n",
             startsWith(as.character(df2plot.1$Line), "D") ~ "Unexposed Pools\n",
             startsWith(as.character(df2plot.1$Line), "E") ~ "Unexposed Pools\n",
             startsWith(as.character(df2plot.1$Line), "F") ~ "Unexposed Pools\n"
  ))

df2plot.1 <- df2plot.1 %>% 
  mutate(
         X.Labs = if_else(Line == "SG200", "", as.vector(Line)))

df2plot.1$Strip.Labs <- factor(df2plot.1$Strip.Labs, 
                             levels = c("Initial\nStrain", "Exposed Pools\n", "Unexposed Pools\n"))


df2plot.1$H2O2 <- gsub("([0-9]+)(mM)", "\\1 \\2", df2plot.1$H2O2)

plot.Figure.1 <- df2plot.1 %>% 
  ggplot() +
  geom_col(aes(x = X.Labs, y = CFU_Ratio_Mean, fill = H2O2, col = H2O2),
           position = position_dodge(width = 0.92)) +
  #geom_hline(yintercept = -0.10, linewidth = 0.7) +
  scale_y_continuous(labels = percent, limits = c(-0.10,1), breaks = c(0, 0.25, 0.50, 0.75, 1.0)) +
  theme_classic() +
  scale_fill_manual(values = c("gray60", "black"))+ 
  scale_color_manual(values = c("gray60", "black"))+ 
  labs(y = "Percentage of Surviving Cells\n") +
  facet_grid(~Strip.Labs, scales = "free", space = "free", switch = "both")+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = c(0.9, 0.8),
    strip.background = element_blank(),
    # axis lines
    axis.line.x = element_blank(),
    #axis.line.y = element_blank(),
    
    axis.title.x = element_blank(),
    axis.text.x = element_text(color = "black", vjust = 34.5),
    
    axis.ticks.x = element_blank(),
    strip.text = element_text(vjust = 2, size = 11)
    ) +
  guides(fill = guide_legend(title = expression("[H"["2"]*"O"["2"]*"] Shock")),
         color = "none"); plot.Figure.1

## save the plot

dirSavePlots <- "Figures"

if ( dir.exists(dirSavePlots) ){
  print (paste0('Directory: "',  dirSavePlots, '" Already Exists!!'))
} else {
  print (paste0('Directory: "',  dirSavePlots, '" does not Exists!!'))
  print (paste0('Directory: "',  dirSavePlots, '" will be Created!!'))
  dir.create(dirSavePlots)
}

file.Figure.1 <- paste(dirSavePlots, "Figure1_USMA_Paper.png", sep = "/")
getwd()
ggsave(filename = file.Figure.1, plot = plot.Figure.1, width = 13, height = 9, dpi = 300, units = "cm")



rm(list = ls())
