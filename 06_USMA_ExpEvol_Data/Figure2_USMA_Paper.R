# author: jcuamatzi

# load libraries
libraries <- c("data.table", "dplyr", "ggplot2", "cowplot", "rstatix", "ggpubr", "tidyr", "svglite")

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

# read data set

df <- fread("F02_ExpEvol_Inheritance_Data/Halos_ResistanceInheritance_withoutH2O2.csv")

df$Lineages <- gsub("L", "Line ", df$Lineages)


#### PLOTTING
df$Group <- factor(df$Group, levels = c("SG200", "H2O2_exposed", "Control"))

df <- df %>% 
  mutate(Strip.Labs = if_else(Lineages == "Initial", " ", as.vector(Lineages)))

#head(df)

df <- df %>% mutate(Fill = 
                case_when(Day == "1 day" & Strain == "SG200" ~ "SG200",
                          Day == "7 day" & Strain == "SG200" ~ "SG200 + 35 Gen.",
                          Day == "1 day" & Strain != "SG200" ~ "200 Generations",
                          Day == "7 day" & Strain != "SG200" ~ "235 Generations"))


#levels(factor(df$Fill))

df$Fill <- factor(df$Fill, levels = c("SG200", "SG200 + 35 Gen.",
                                      "200 Generations", "235 Generations"))

#df
#df <- df %>% mutate(Strip.Labs.2 =  if_else(Strip.Labs %in% c("Line B", "Line E"), "Exposed", "")) 

df <- df %>% mutate(Strip.Labs.2 = case_when(
  startsWith(df$Lineages, "Initial") ~ "Initial\nStrain",
  startsWith(df$Lineages, "Line A") ~ "",
  startsWith(df$Lineages, "Line B") ~ "Exposed Colonies\n",
  startsWith(df$Lineages, "Line C") ~ "",
  startsWith(df$Lineages, "Line D") ~ "",
  startsWith(df$Lineages, "Line E") ~ "Unexposed Colonies\n",
  startsWith(df$Lineages, "Line F") ~ ""
))

df %>% dplyr::filter(Strain != "SG200") %>% dplyr::group_by(Group) %>% t_test(cm ~ Day, paired = T)

#plot.Figure.2 <- df %>% group_by(Day, Strain, Lineages, Group, Strip.Labs, Fill) %>% 
plot.Figure.2 <- 
  df %>% group_by(Day, Strain, Lineages, Group, Strip.Labs, Strip.Labs.2, Fill) %>% 
  summarise(Mean = mean(cm)) %>% 
  mutate(InhibitionArea = (pi*(Mean/2)^2)*10) %>%  # express the inhibition area in mm^2
  ungroup() %>% 
  ggplot()+
  geom_col(aes(x = Strain, y = InhibitionArea,
               fill = Fill, linetype = Fill, color = Fill), 
           position = position_dodge(width = 0.80), 
           alpha = 0.8, width = 0.80, linewidth = 0.30) +
  
  scale_fill_manual(values = c("gray90", "gray90", "gray60", "black")) +
  
  scale_color_manual(values = c("black", "black","gray60", "black")) +
  
  scale_linetype_manual(values = c("solid", "dotted", "blank", "blank"))+
  
  labs(
    #x = "Colonies", 
       y = expression("Zone of Inhibition by H"["2"]*"O"["2"]*" (mm"^"2"*")")) +
  #facet_grid(~ as.factor(Strip.Labs) + as.factor(Strip.Labs.2 ),
  facet_grid(~ as.factor(Strip.Labs) + Strip.Labs.2,
             scales = "free", space = "free", switch = "both")+
  theme_classic()+
  scale_x_discrete(labels = c("SG200" = "",
                              "T20.L1.C1" = "1", "T20.L1.C2" = "2",
                              "T20.L1.C3" = "3", "T20.L1.C4" = "4", 
                              "T20.L1.C5" = "5",
                              "T20.L2.C1" = "1", "T20.L2.C2" = "2",
                              "T20.L2.C3" = "3", "T20.L2.C4" = "4", 
                              "T20.L2.C5" = "5",
                              "T20.L3.C1" = "1", "T20.L3.C2" = "2",
                              "T20.L3.C3" = "3", "T20.L3.C4" = "4", 
                              "T20.L3.C5" = "5",
                              
                              "U20.L1.C1" = "1", "U20.L1.C2" = "2",
                              "U20.L2.C1" = "1", "U20.L2.C2" = "2",
                              "U20.L3.C1" = "1", "U20.L3.C2" = "2" )) +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0, 120, 20)) +
  
  theme(
    # legend aes
    legend.position = c(0.5, .9),
    legend.direction = "horizontal",
    legend.title = element_text(colour="white"),
    # 
    axis.line.x = element_blank(),
    axis.text.x = element_text(vjust = 6, color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    
    # strip
    strip.background = element_blank(),
    
    strip.text = element_text(size = 11, color = "black"),
    strip.placement = "outside",
    
    # axis titles
    axis.title.y = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(),
    
    ) +
  guides(fill = guide_legend(title = "Generations"),
         color = guide_legend(title = "Generations"),
         linetype = guide_legend(title = "Generations")
         ); plot.Figure.2

## save the plot

dirSavePlots <- "Figures"

if ( dir.exists(dirSavePlots) ){
  print (paste0('Directory: "',  dirSavePlots, '" Already Exists!!'))
} else {
  print (paste0('Directory: "',  dirSavePlots, '" does not Exists!!'))
  print (paste0('Directory: "',  dirSavePlots, '" will be Created!!'))
  dir.create(dirSavePlots)
}

file.Figure.2 <- paste(dirSavePlots, "Figure2_USMA_Paper.New.svg", sep = "/")

ggsave(filename = file.Figure.2, plot = plot.Figure.2,
       width = 25, height = 12, dpi = 300, units = "cm")


rm(list = ls())

