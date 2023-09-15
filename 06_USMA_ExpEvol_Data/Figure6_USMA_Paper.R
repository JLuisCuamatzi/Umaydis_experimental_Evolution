# author: jcuamatzi
 # code to plot the infection symptoms was adapted from https://gitlab.gwdg.de/molsysevol/umag_11064/-/blob/master/Infections/InfectionResultsAnalysis.R?ref_type=heads
# load libraries
libraries <- c("ggplot2", "dplyr", "tidyr", "scales", "data.table", "RColorBrewer", "ggthemes", "cowplot", "wesanderson", "rstatix")

for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    suppressPackageStartupMessages(install.packages(lib, dependencies = TRUE))
  }
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

# set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# clean the env
rm(list = ls())

# Define color palette for SG200, T20.LB1, and T20.LC.1
color.sg <- "black"
color.lb <- "#7F6000"
color.lc <- "#C55A11"

## plot oexcat
# Read file with data of CFU in oexUMAG_11067
df.oex <- fread("oexUMAG_11067/UMAG_11067_oex_Mutant.csv")

df.oex.mean <- df.oex %>% 
  group_by(Strain, H2O2) %>% 
  summarise(MeanCFU = mean(CFU)) %>% 
  ungroup()

df.oex.mean <- df.oex.mean %>% 
  mutate(CFU_Type = case_when(
    startsWith("0 mM", df.oex.mean$H2O2) ~ "CFU_Reference",
    startsWith("10 mM",df.oex.mean$H2O2) ~ "CFU_Target"))


df.oex.mean$H2O2 <- gsub(" ", "", df.oex.mean$H2O2)

# Pivot the table
df.oex.mean <- df.oex.mean %>% group_by(Strain) %>% 
  reframe(CFU_0mM = MeanCFU[H2O2 == "0mM"],
            CFU_10mM = MeanCFU[H2O2 == "10mM"]) %>%
  ungroup()

# estimate the percentage of cfu
# we calculate 10 mM against 0 mM in each group
df.oex.mean <- df.oex.mean %>%
  mutate(across(c(`CFU_0mM`:`CFU_10mM`),
                .fns = ~./`CFU_0mM`,
                .names = "{.col}_Ratio"))

# Pivot from wider to longer
df.oex.mean <- df.oex.mean %>%
  select(Strain, CFU_0mM_Ratio: CFU_10mM_Ratio) %>%
  pivot_longer(cols = CFU_0mM_Ratio: CFU_10mM_Ratio, names_to = "H2O2", values_to = "CFU_Ratio") %>%
  setDT()

df.oex.mean$H2O2 <- gsub("CFU_", "", df.oex.mean$H2O2)    # Remove 'CFU_'
df.oex.mean$H2O2 <- gsub("_Ratio", "", df.oex.mean$H2O2)  # Remove '_Ratio'

df.oex.mean$Strain <- factor(df.oex.mean$Strain, levels = c("SG200", "oexCAT","T20.LC.1"))

df.oex <- df.oex %>% group_by(Strain) %>%  arrange(H2O2) %>% ungroup()

df.oex <- df.oex %>% mutate(CFU_Type = case_when(
  startsWith(df.oex$H2O2, "0 mM") ~ "CFU_Reference",
  startsWith(df.oex$H2O2, "10 mM") ~ "CFU_Target"
))

df.oex <- df.oex %>% group_by(Strain, Replicate, H2O2) %>% summarise(MeanCFU = mean(CFU) ) %>% 
  ungroup()

df.oex <- df.oex %>% mutate(CFU_Type = case_when(
  startsWith(df.oex$H2O2, "0 mM") ~ "CFU_Reference",
  startsWith(df.oex$H2O2, "10 mM") ~ "CFU_Target"
))

df.oex$Percentage <- 0

for (i in 1:length(df.oex$Strain)){
  if (df.oex[i, 5] == "CFU_Reference"){
    df.oex[i, 6] <- df.oex[i, 4] / df.oex[i, 4]
    df.oex[i + 1, 6] <- df.oex[i + 1, 4] / df.oex[i, 4]
    
  }
}

df.oex <- df.oex %>% filter(H2O2 == "10 mM")
df.oex$Strain <- factor(df.oex$Strain, levels = c("SG200", "oexCAT", "T20.LC.1"))

df.oex.mean$H2O2 <- gsub("([0-9]+)(mM)", "\\1 \\2", df.oex.mean$H2O2)
df.oex.mean <- df.oex.mean[H2O2 != "0 mM"]
df.oex.mean$Strain <- factor(df.oex.mean$Strain, levels = c("SG200", "oexCAT", "T20.LC.1"))

# Figure 6.A
plot.oex.mutant <- ggplot() +
  geom_point(data = df.oex, aes(x = Strain, y = Percentage, color = Strain), 
             size = 3, alpha = 0.95) +
  geom_col(data = df.oex.mean, aes(x = Strain, y = CFU_Ratio, fill = Strain), color = "black", alpha = 0.75)+
  scale_y_continuous(labels = percent, limits = c(0,1)) +
  scale_x_discrete(labels = c("SG200" = "Initial \nStrain",
                              "oexCAT" = "oex \nUMAG_11067",
                              "T20.LC.1" = "T20. \nLC.1"))+
  scale_fill_manual(values = wes_palette("Cavalcanti1")) +
  #scale_fill_colorblind() +
  #scale_fill_manual(values = c("black", "darkgreen", color.lc))+
  scale_color_manual(values = wes_palette("Cavalcanti1"))+
  #scale_color_colorblind() +
  
  labs(y = "Percentage of Surviving Cells\n", x = "\nStrain")+
  theme_classic() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 12, color = "black"),
        #axis.text.x = element_blank(),
        #strip.text = element_text(size = 14, face = "bold", color = "black"),
        #
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"));plot.oex.mutant


# read the data frame with information about USMA infection on maize plants
df <- fread("oexUMAG_11067/oexCAT_Maize_Infection_11dpi.csv") # fread from data.table (data.table::fread())
df <- df[, -12]

# create a data frame in which each symptom has a raking from 1 to the highest symptom
symptom.weight.df <- data.frame("Symptom" = names(df[, 4:ncol(df)]),
                                "Weight" = seq(1, length(4:ncol(df)), 1))

# create a column with a unique id for each plant (paste strain + plantNum)
df$Plant_ID <- paste(df$Strain, df$`#Plant`, sep = "_")

# remove columns without relevant information (#Pot, #Plant, #PlantStrain)
df <- df[ ,-c(1:2) ]

## Transform the data frame to create a column called 'Symptom_name' in which is listed all symptoms 
# also, the values 1 = present, 0 = absent will be in a column called 'Presence'
# we used pivot_longer() from tidyr (tidyr::pivot_longer())

df.2 <- df %>% pivot_longer(names_to = "Symptom_name", # name for column with each symptom
                            values_to = "Presence", # name for column with information about presence/absence
                            cols = 2:9) %>%  # colums with information about presence/absence 
  left_join(select(symptom.weight.df, Symptom, Weight), # add the ranking for each symptom
            by = c("Symptom_name" = "Symptom")) 

## For the analysis, we use the highest symptom in the plant, so, we need to identify the highest symptom

# 
df.2 <- df.2 %>% 
  # here, we multiply the weight of each symptom by the value in presence (if a symptom is absent, then the result will be 0) 
  mutate(Weight_Symptom = df.2$Presence * df.2$Weight )%>% 
  group_by(Strain, Plant_ID) %>% # group by Strain, and Plant_ID to know the highest symptom in each plant
  summarise (MaxSymptom = max(Weight_Symptom, na.rm = T)) %>% # to summarise by plant (we keep only the max value)
  ungroup() %>% # ungroup for next functions! 
  left_join(select(symptom.weight.df, Symptom, Weight), # add the name of each symptom
            by = c("MaxSymptom" = "Weight"))

# count the frequency of each symptom

df.3 <- df.2 %>% group_by(Strain, Symptom) %>% 
  count(Symptom) %>% 
  pivot_wider(names_from = Symptom,
              values_from = n) %>% 
  ungroup() %>% 
  mutate_all(~replace(., is.na(.),0)) 

# df.3[symptom.weight.df$Symptom[!(symptom.weight.df$Symptom %in% colnames(df.3))]] <- 0
# df.3

df.3.names <- df.3[,1]

df.3.symp <- df.3[,-1]
df.3.symp <- df.3.symp[, order(colnames(df.3.symp))]

df.3 <- bind_cols(df.3.names, df.3.symp)

rm(df.3.names, df.3.symp)

#end_time <- Sys.time()
#end_time - start_time

df.3
# chi-sq test for SG200 vs oexCAT
df.3.SG200vsxCAT <- as.matrix(df.3[c(5,4), -1])
df.3.SG200vsxCAT
#p.sg200.oexcat <- fisher.test(df.3.SG200vsxCAT, simulate.p.value = TRUE, B = 1000000)  ## p-value = 0.00036
chisq.test(df.3.SG200vsxCAT, simulate.p.value = TRUE, B = 1000000)   ## p-value = 0.00029
#p.sg200.oexcat$p.value

# chi-sq test for SG200 vs LC.1
df.3.SG200vsLC.1 <- as.matrix(df.3[c(4,2), -1])

chisq.test(df.3.SG200vsLC.1, simulate.p.value = TRUE, B = 1000000)  ## p-value = 0.005038
#p.sg200.lc.1$p.value

# chi-sq test for oexCAT vs LC.1
df.3.xCATvsLC.1 <- as.matrix(df.3[c(5,2), -1])
df.3.xCATvsLC.1 

#chisq.test(df.3.xCATvsLC.1, simulate.p.value = TRUE, B = 1000000)  ## here I removed Chlorosis (0
chisq.test(df.3.xCATvsLC.1[,c(1,3:8)], simulate.p.value = TRUE, B = 1000000)  ## here I removed Chlorosis (0 observations in both data sets). P = 0.2993

rm(df.2, symptom.weight.df)

# proportion of each symptom 
df.4 <- df.3 %>% 
  tidyr::pivot_longer(names_to = "Symptom",
                      values_to = "SymptomFreq",
                      cols = 2:ncol(df.3))%>% 
  dplyr::group_by(Strain) %>%
  dplyr::mutate(Percentage = SymptomFreq/sum(SymptomFreq)) %>% setDT()

#
df.3.tmp <- df.3[c(5, 4, 2) ,]

counts <- data.frame(Strain = c("SG200","oexCAT", "LC.1"),
                     Count = paste("n =", rowSums(df.3.tmp[,-1])))

# colors for plot
colors <- brewer.pal(7, "YlOrBr")

df.4$Strain <- factor(df.4$Strain, levels = c("SG200", "oexCAT", "LC.1"))

# no H2O and no inoculo

df.4 <- df.4[Strain != "H2O"]
df.4 <- df.4[Strain != "NoInoculo"]

df.4$Symptom <- gsub("[1-9].", "", gsub("_", " ", df.4$Symptom))

df.4$Symptom <- gsub("No symptom", "No Symptoms", df.4$Symptom)

df.4$Symptom <- factor(df.4$Symptom, levels = c("No Symptoms", "Chlorosis", "Ligula swelling",
                                                "Small tumors", "Normal tumors", "Heavy tumors stem", 
                                                "Heavy tumors", "Dead Plant"))
levels(factor(df.4$Strain))

# Figure 6.B
plot.infection <- ggplot() + 
  geom_col(data = df.4, 
           aes(x = Strain, 
               y = Percentage, 
               fill = Symptom), alpha = 0.9, col = "black", width = 0.6) +
  scale_x_discrete(labels = c("SG200" = "Initial \nStrain", 
                              "oexCAT" = "oex \nUMAG_11067", 
                              "LC.1" = "T20. \nLC.1")) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("No Symptoms" = "white",
                               "Chlorosis" = colors[1],
                               "Ligula swelling" = colors[2],
                               "Small tumors" = colors[3],
                               "Normal tumors" = colors[4],
                               "Heavy tumors stem" = colors[5],
                               "Heavy tumors" = colors[6],
                               "Dead Plant" = colors[7])) +
  geom_text(data = counts, aes(x = Strain, y = 1.05, label = Count), size = 3) +
  labs(y = "Percentage of Infected Plants", x = "\nStrain")+
  theme_classic()+ 
  theme(
    # legend
    #legend.position = c(,0.7),
    legend.text = element_text(size = 8),
    #legend.key.size = unit(0.1, "cm"),
    #legend.position = "top",
    
    # axis
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text( size = 10, color = "black"),
        #strip
        strip.background = element_blank(),
        strip.text = element_blank());plot.infection


Figure.6 <- plot_grid(plot.oex.mutant, plot.infection, 
                      rel_widths = c(0.75, 1), scale = 0.9, labels = c("A)", "B)"))

Figure.6

#### export plots
dirSavePlots <- "Figures"

if ( dir.exists(dirSavePlots) ){
  print ("Directory already exists!!")
} else {
  dir.create(dirSavePlots)
}

plot.Figure6 <- paste(dirSavePlots, "Figure6_USMA_Paper.png", sep = "/")

ggsave(filename = plot.Figure6, plot = Figure.6,
       width = 22, height = 10, units = "cm", dpi = 300)


rm(list = ls())
