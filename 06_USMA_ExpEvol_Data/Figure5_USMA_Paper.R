# created by j cuamatzi
# script for figure 5

# load libraries
libraries <- c("data.table", "readxl","dplyr", "ggplot2", "tidyr", "ggpubr", "cowplot", "scales", "rstatix")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}


# clean env
rm(list = ls())


# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## Figure 5.A
df <- fread("F05_ExpEvol_Resistance_Chr9_Data/InhibitionEachColony.csv")

df <- df %>% 
  mutate(X.Labs = case_when(
    endsWith(df$Strain, "SG200") ~ "SG200",
    endsWith(df$Strain, ".1") ~ "1",
    endsWith(df$Strain, ".2") ~ "2",
    endsWith(df$Strain, ".3") ~ "3",
    endsWith(df$Strain, ".4") ~ "4",
    endsWith(df$Strain, ".5") ~ "5"))

df.col <- df[HaloSize_cm != 0] # remove control plates (cm = 0)

# indicate levels in strain for plot
df.col$Strain <- factor(df.col$Strain, 
                        levels = c("SG200", 
                                   "U20.LD.1", "U20.LD.2", "U20.LE.1", "U20.LE.2", "U20.LF.1", "U20.LF.2",
                                   "T20.LA.1","T20.LA.2","T20.LA.3", "T20.LA.4","T20.LA.5",
                                   "T20.LB.1","T20.LB.2","T20.LB.3", "T20.LB.4","T20.LB.5",
                                   "T20.LC.1","T20.LC.2","T20.LC.3", "T20.LC.4","T20.LC.5"))

# indicate levels in lines for plot
df.col$Line <- factor(df.col$Line, levels = c("Initial", "D", "E", "F", "A", "B", "C"))



# SG200
SG200.Halo <- df.col[Line == "Initial"]

SG200.Halo$ZoneInhibition <- (pi*((SG200.Halo$HaloSize_cm/2)^2)*10)

mean.SG200 <- mean(SG200.Halo$ZoneInhibition)
sd.SG200 <- sd(SG200.Halo$ZoneInhibition)


df.col.2 <- df.col[!Line %in% c("Initial", "D", "E", "F")]
df.col.2 <- df.col.2 %>% mutate(Chr9State = 
                                  case_when(
                                    endsWith(df.col.2$Chr9.LA, "1X") ~ "Euploidy",
                                    endsWith(df.col.2$Chr9.LA, "2X") ~ "Aneuploidy",
                                    endsWith(df.col.2$Chr9.LA, "3X") ~ "Aneuploidy"
                                  )) 

df.col.2$Chr9State <- factor(df.col.2$Chr9State, levels = c("Euploidy", "Aneuploidy"))
df.col.2
df.col.2.stats <- df.col.2 %>% group_by(Line) %>%  t_test(HaloSize_cm ~ Chr9State) %>% 
  adjust_pvalue() %>% add_significance()


df.col.2.stats.2 <- df.col.2 %>% t_test(HaloSize_cm ~ Chr9.LA) %>% adjust_pvalue() %>% add_significance()

df.col.2.stats.2

my_comparisons <- list( c("1X", "2X"), c("1X", "3X"), c("2X", "3X") )

df2plot <- df.col.2 %>% mutate(InhibitionArea = ((pi)*(HaloSize_cm/2)^2)) %>% 
  group_by(Chr9.LA, Line, Strain) %>% 
  summarise(HaloAreaMean = mean(InhibitionArea),
            HaloAreaSD = sd(InhibitionArea)) %>% 
  ungroup()

df2plot$Strain <-  (df2plot$Strain)

df2plot$Line <- paste("Line ", df2plot$Line, sep = "")

df2plot$Chr9.LA <- factor(df2plot$Chr9.LA, levels = c("1X", "2X", "3X"))

# add dots to each bar

df2plot.dots <- df.col.2 %>% mutate(InhibitionArea = ((pi)*(HaloSize_cm/2)^2))

df2plot.dots$Line <- paste("Line ", df2plot.dots$Line, sep = "")

df2plot <- df2plot %>% arrange(Strain) %>% mutate(X.Labs = rep(seq(1,5,1), 3) ) 


 
Figure.5A <- ggplot(data = df2plot) +
  geom_col( aes(x = X.Labs, y = (HaloAreaMean*10), fill = Chr9.LA, color = Chr9.LA), 
           alpha = 0.75, linewidth = 0.15)+
  geom_point(data = df2plot.dots, aes(x = as.numeric(X.Labs), y = (InhibitionArea*10)), color = "black", alpha = 0.5)+
  facet_grid(~Line, scales = "free", space = "free", switch = "both") +
  geom_hline(yintercept = 0)+
  theme_classic() +
  #scale_y_continuous(limits = c(0, (6.5*10)), breaks = seq(0, (6*10), (1*10))) +
  scale_y_continuous(limits = c(0, (12*10)), breaks = seq(0, (12*10), (2*10))) +
  labs(x ="Colonies at Generation 200", 
       y = expression("Inhibition Zone by H"["2"]*"O"["2"]*" (mm"^"2"*")")) +
  scale_fill_manual(values = c( 
    "white", #1X
    "gray", # 2X
    "black"  # 3X
    )) +
  scale_color_manual(values = c( 
    "black", #1X
    "black", # 2X
    "black"  # 3X
  )) +  
  theme(
    legend.position = c(0.5, 0.94),
    legend.direction = "horizontal",
    legend.key.size = unit(0.35, "cm"),
    
    
    axis.title.x = element_text(margin = unit(c(3,0,0,0), "mm")),
    axis.title.y = element_text(margin = unit(c(0,5,0,0), "mm")),
    
    axis.line.x = element_blank(),
    axis.text.x = element_text(vjust = 15, color = "black", size = 9),
    axis.ticks.x = element_blank(),
    #
    axis.text.y = element_text(color = "black", size = 9),
    #
    strip.background = element_blank(),
    strip.text = element_text(size = 11, color = "black"),
    #
    axis.title = element_text(size = 12, color = "black")) +
  guides(fill = guide_legend(title = "Copy Number\nLeft Arm Chr 9",
                             title.theme = element_text(size = 8, face = "bold")),
         color = guide_legend(title = "Copy Number\nLeft Arm Chr 9",
                              title.theme = element_text(size = 7, face = "bold")))+
  geom_hline(aes(yintercept = mean.SG200), color = "red", linetype = "dotted" );Figure.5A




  ## Figure 5.B
df.qpcr <- read_excel(path = "F05_ExpEvol_Resistance_Chr9_Data/UMAG_11067_qPCR.xlsx", sheet = "Fold_Change") # read the file
df.qpcr <- setDT(df.qpcr) # transform to data table

stat.test <- df.qpcr %>% 
  group_by(Condition) %>% 
  t_test(FoldChange ~ Strain) %>% 
  add_xy_position(x = "Strain")

stat.test <- stat.test[stat.test$Condition != "0 mM",]

stat.test$label <- sprintf("%.3f", stat.test$p)

# Estimate the mean of the fold change

df.mean <- df.qpcr %>% group_by(Strain, Condition) %>% 
  summarise(Mean.FC = mean(FoldChange)) %>% 
  ungroup() %>% setDT()

df.mean <- df.mean[Condition != "0 mM"]
df.qpcr <- df.qpcr[Condition != "0 mM"]



df.mean <- df.mean %>% mutate(Chr9.LA = case_when(
  startsWith(df.mean$Strain, "SG200") ~ "1X",
  startsWith(df.mean$Strain, "T20.LB.1") ~ "2X",
  startsWith(df.mean$Strain, "T20.LC.1") ~ "3X"
))


df.qpcr <- df.qpcr %>% mutate(Chr9.LA = case_when(
  startsWith(df.qpcr$Strain, "SG200") ~ "1X",
  startsWith(df.qpcr$Strain, "T20.LB.1") ~ "2X",
  startsWith(df.qpcr$Strain, "T20.LC.1") ~ "3X"
))

# box plot con ggplot
Figure.5B <- ggplot()+
  #geom_boxplot()+
  geom_col(data = df.mean, 
           aes(x = Strain, y = Mean.FC, 
               fill = Chr9.LA, color = Chr9.LA),
           alpha = 0.75, linewidth = 0.15) +
  geom_point(data = df.qpcr, aes(x = Strain, y = FoldChange), 
             size = 2, alpha = 0.5, color = "black")+
  geom_hline(yintercept = 0)+
  facet_grid(~Condition) +
  theme_classic()+
  scale_y_continuous(limits = c(0, 5), breaks = seq(0,5, 1))+
  scale_x_discrete(labels = c("SG200" = "Initial \nStrain", 
                              "T20.LB.1" = "1 \nLine B",
                              "T20.LC.1" = "1 \nLine C"))+
  labs(x = "Strain", y = "Fold Change in Catalase Gene Expression")+
  ## colors
  scale_fill_manual(values = c( 
    "white", #1X
    "gray", # 2X
    "black"  # 3X
  )) +
  scale_color_manual(values = c( 
    "black", #1X
    "black", # 2X
    "black"  # 3X
  )) +  
  
  ##
  theme(
    # legend
    legend.position = "none",
    
    # axis title
        
        axis.title.x = element_text(margin = unit(c(3,0,0,0), "mm")),
        axis.title.y = element_text(margin = unit(c(0,5,0,0), "mm")),
        
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = 5, color = "black", size = 9),
        
        axis.title = element_text(size = 12, color = "black"),
        #axis.text.x = element_blank(),
        #strip.text = element_text(size = 14, face = "bold", color = "black"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        
        axis.text.y = element_text(size = 9, color = "black")) +
  stat_pvalue_manual(data = stat.test, label = "label")


Figure.5 <- plot_grid(Figure.5A, Figure.5B, rel_widths = c(1.0, 0.5), labels = c("A)", "B)"), scale = 0.9)

#Figure.5

dirSavePlots <- "Figures"

if ( dir.exists(dirSavePlots) ){
  print ("Directory already exists!!")
} else {
  dir.create(dirSavePlots)
}

plot.Figure5 <- paste(dirSavePlots, "Figure5_USMA_Paper.svg", sep = "/")


ggsave(filename = plot.Figure5, plot = Figure.5,
       width = 18, height = 12, units = "cm", dpi = 300)



rm(list = ls())


