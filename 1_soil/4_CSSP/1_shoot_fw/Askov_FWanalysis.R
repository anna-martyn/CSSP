options(warn=-1)

# cleanup
rm(list=ls())

# directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# files
weight.file <- paste("Input_OX&SYM_mutant_harvest_2.txt", sep = "")

# load data
weight <- read.table(weight.file, header=T, sep="\t")
weight <- as.data.frame(weight)
weight$Fresh_weight <- gsub("\\,",".",weight$Fresh_weight)
weight$Fresh_weight <- as.numeric(weight$Fresh_weight)

#Filter data to remove OX info
idx <- weight$Genotype%in% c("WT","symrk","ccamk","nsp1","nsp2")
weight_filt <- weight[idx,]

idx <- weight_filt$Experiment%in% c("Lotus_CSSP","Barley_CSSP")
weight_filt_CSSP <- weight_filt[idx,]

idx <- weight_filt_CSSP$Plant_species%in% c("Barley")
Barley_weight <- weight_filt_CSSP[idx,]

idx <- weight_filt_CSSP$Plant_species%in% c("Lotus")
Lotus_weight <- weight_filt_CSSP[idx,]

# set color for the genotype
# colors <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
#                      color=c("#33a02c","#ff7f00","#1f78b4","#e31a1c","#ffd700"))

colors <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
                     color=c("#A9C289","#FEDA8B","#FDB366","#C0E4EF","#6EA6CD"))

# colors <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
#                      color=c("#A9C289","#FEDA8B","#FDB366","#98CAE1","#4A7BB7"))


# set order for plants
Lotus_weight$Soil_type <- factor(Lotus_weight$Soil_type, levels = c("NPK","PK","UF"))
Lotus_weight$Genotype <- factor(Lotus_weight$Genotype, levels = colors$group)

Barley_weight$Soil_type <- factor(Barley_weight$Soil_type, levels = c("NPK","PK","UF"))
Barley_weight$Genotype <- factor(Barley_weight$Genotype, levels = colors$group)

# load plotting functions
library("ggplot2")
library("scales")
library("grid")
library(RColorBrewer)

## boxplot
dodge <- position_dodge(width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1),
                    axis.text.y = element_text(size = 20),
                    #legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

p1 <- ggplot(Lotus_weight, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  #geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  geom_boxplot(aes(shape=Soil_type), width=0.5,position = dodge,outlier.color = NA)+
  geom_jitter(aes(shape=Soil_type), position=position_jitterdodge(jitter.width = 0.3), size=2, alpha=0.3)+
  scale_fill_manual(values=as.character(colors$color)) +
  facet_wrap(~Soil_type, scales = "free_x", nrow=1)+
  main_theme +
  ylab("Shoot fresh weight/plant [g]")+
  scale_y_continuous(limits = c(0, 0.15))+
  ggtitle("Lotus")+
  theme(legend.position="none", 
        plot.title = element_text(size = 20, face="bold"),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1, colour = "black"),
        legend.key.size = unit(1,"cm"))
p1

ggsave(paste("Lotus_CSSP_fw_BOX_v2.png", sep=""), p1, width=6, height=5)
ggsave(paste("Lotus_CSSP_fw_BOX_v2.pdf", sep=""), p1, width=6, height=5)

p2 <- ggplot(Barley_weight, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  #geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  geom_boxplot(aes(shape=Soil_type), width=0.5,position = dodge,outlier.color = NA)+
  geom_jitter(aes(shape=Soil_type), position=position_jitterdodge(jitter.width = 0.3), size=2, alpha=0.3)+
  scale_fill_manual(values=as.character(colors$color)) +
  facet_wrap(~Soil_type, scales = "free_x", nrow=1)+
  main_theme +
  ylab("Shoot fresh weight/plant [g]")+
  scale_y_continuous()+
  ggtitle("Barley")+
  theme(legend.position="none", 
        plot.title = element_text(size = 20, face="bold"),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1, colour = "black"),
        legend.key.size = unit(1,"cm"))
p2

ggsave(paste("Barley_CSSP_fw_BOX_v2.png", sep=""), p2, width=6, height=5)
ggsave(paste("Barley_CSSP_fw_BOX_v2.pdf", sep=""), p2, width=6, height=5)

#Subset by soil type for each plant species
idx <- Lotus_weight$Soil_type%in% c("NPK")
Lotus_weight_NPK <- Lotus_weight[idx,]

idx <- Lotus_weight$Soil_type%in% c("PK")
Lotus_weight_PK <- Lotus_weight[idx,]

idx <- Lotus_weight$Soil_type%in% c("UF")
Lotus_weight_UF <- Lotus_weight[idx,]

idx <- Barley_weight$Soil_type%in% c("NPK")
Barley_weight_NPK <- Barley_weight[idx,]

idx <- Barley_weight$Soil_type%in% c("PK")
Barley_weight_PK <- Barley_weight[idx,]

idx <- Barley_weight$Soil_type%in% c("UF")
Barley_weight_UF <- Barley_weight[idx,]

# statistical analysis (separate ANOVA and TukeyHSD for each Plant-Soil combination)
library(car)
library(multcompView)

####LOTUS
## NPK
ano <- aov(Fresh_weight ~ Genotype, data=Lotus_weight_NPK)
anova(ano)
pairwise <- TukeyHSD(ano)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}

pairwise[["Genotype"]] <- na.omit(pairwise[["Genotype"]])
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file =  "Lotus_CSSP_NPK_fw_ANOVA.csv")

## PK
ano <- aov(Fresh_weight ~ Genotype, data=Lotus_weight_PK)
anova(ano)
pairwise <- TukeyHSD(ano)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}

pairwise[["Genotype"]] <- na.omit(pairwise[["Genotype"]])
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file =  "Lotus_CSSP_PK_fw_ANOVA.csv")

## UF
ano <- aov(Fresh_weight ~ Genotype, data=Lotus_weight_UF)
anova(ano)
pairwise <- TukeyHSD(ano)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}

pairwise[["Genotype"]] <- na.omit(pairwise[["Genotype"]])
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file =  "Lotus_CSSP_UF_fw_ANOVA.csv")

####BARLEY
## NPK
ano <- aov(Fresh_weight ~ Genotype, data=Barley_weight_NPK)
anova(ano) 
pairwise <- TukeyHSD(ano)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}

pairwise[["Genotype"]] <- na.omit(pairwise[["Genotype"]])
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file =  "Barley_CSSP_NPK_fw_ANOVA.csv")

## PK
ano <- aov(Fresh_weight ~ Genotype, data=Barley_weight_PK)
anova(ano) 
pairwise <- TukeyHSD(ano)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}

pairwise[["Genotype"]] <- na.omit(pairwise[["Genotype"]])
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file =  "Barley_CSSP_PK_fw_ANOVA.csv")

## UF
ano <- aov(Fresh_weight ~ Genotype, data=Barley_weight_UF)
anova(ano)
pairwise <- TukeyHSD(ano)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}

pairwise[["Genotype"]] <- na.omit(pairwise[["Genotype"]])
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file =  "Barley_CSSP_UF_fw_ANOVA.csv")

