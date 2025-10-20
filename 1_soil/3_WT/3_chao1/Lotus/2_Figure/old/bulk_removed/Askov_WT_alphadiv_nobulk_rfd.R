options(warn=-1)

# cleanup
rm(list=ls())

setwd("O:/Nat_MBG-PMg/Anna_Martyn/0_MARTYN_THORSGAARD/2_Soil_experiment/2_WT/Anna analysis separate/3_chao1/Lotus/rarefied/bulk_removed/")

results.dir <- "O:/Nat_MBG-PMg/Anna_Martyn/0_MARTYN_THORSGAARD/2_Soil_experiment/2_WT/Anna analysis separate/3_chao1/Lotus/rarefied/bulk_removed/"
figures.dir <- "O:/Nat_MBG-PMg/Anna_Martyn/0_MARTYN_THORSGAARD/2_Soil_experiment/2_WT/Anna analysis separate/3_chao1/Lotus/rarefied/bulk_removed/"


alpha.file <- paste(results.dir, "Askov_Lotus_10_4_rfd_nobulk_chao1.txt", sep="")
design.file <- paste(results.dir, "Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
alpha <- read.table(alpha.file, sep="\t", header=T, check.names=F, row.names = 1)

### alpha diversity
#shapes <- data.frame(group=c("rhizosphere","root","bulk"), shape=c(16,1,17))
colors <- data.frame(group=c("NPK","PK","UF"), color=c("#1b9e77","#d95f02","#666666"))


# Chao1 index
index <- cbind(alpha[,1], design[match(row.names(alpha), design$SampleID), ])
colnames(index)[1] <- "value"


# boxplots by Compartment

idx <- index$Genotype %in% c( "WT")
index_subset <- index[idx, ]

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggh4x)

index_subset <- index_subset[complete.cases(index_subset),]

l1 <- c("NPK", "PK", "UF")
index_subset$Soil <- factor(index_subset$Soil, levels=l1)
colors <- colors[match(l1, colors$group), ]

#l2 <- c("rhizosphere","endosphere","bulk")
#index_subset$Compartment <- factor(index_subset$Compartment, levels = l2)
#shapes <- shapes[match(l2, shapes$group),]

index_subset$Compartment <- factor(index_subset$Compartment, levels = c("rhizosphere","root","nodules"))


main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))


p15 <- ggplot(index_subset, aes(x=Soil, y=value, fill=Soil)) +
  geom_boxplot(alpha=0.7, position= position_dodge(width = 0.7), outlier.color =NA, width=0.3) +
  geom_jitter(position=position_jitter(0.17), size=4, alpha=1) +
  #geom_jitter(aes(shape=shapes), position=position_jitter(0.17), size=4, alpha=1) +
  scale_fill_manual(values=as.character(colors$color)) +
  facet_nested(~Compartment,scales ="free_x", labeller = label_parsed)+
  #scale_shape_manual(values=shapes$shape) +
  ggtitle("Lotus")+
  labs(x="", y="Chao1 index") +
  main_theme +
  theme(legend.position="none", 
        plot.title = element_text(size = 20, face = "bold"),  # Adjust the size here
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 20),
        legend.key.size = unit(1,"cm"))
p15

ggsave(paste(figures.dir, "Chao1_WT_Lotus_rfd_bulk_removed.png", sep=""), p15, width=7, height=6)


# compute summary statistics

library(dplyr)

idx <- index_subset$Compartment %in% c( "rhizosphere")
index_rhizo <- index_subset[idx, ]

idx <- index_subset$Compartment %in% c( "root")
index_root <- index_subset[idx, ]

idx <- index_subset$Compartment %in% c( "nodules")
index_nodules <- index_subset[idx, ]


chao_summarize1 <- index_rhizo%>%
  group_by(Soil)%>%
  summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value))

chao_summarize2 <- index_root%>%
  group_by(Soil)%>%
  summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value))

chao_summarize3 <- index_nodules%>%
  group_by(Soil)%>%
  summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value))



# statistical analyses by Mann-Whitney U-test

### check the differences one by one
###RHIZO

idx <- index_rhizo$Soil %in% c("NPK","PK")
rhizo_chao_Barley_NPK_PK <- index_rhizo[idx, ]

wilcox.test(value~Soil, data = rhizo_chao_Barley_NPK_PK)

idx <- index_rhizo$Soil %in% c("PK","UF")
rhizo_chao_Barley_PK_UF <- index_rhizo[idx, ]

wilcox.test(value~Soil, data = rhizo_chao_Barley_PK_UF)

idx <- index_rhizo$Soil %in% c("NPK","UF")
rhizo_chao_Barley_NPK_UF <- index_rhizo[idx, ]

wilcox.test(value~Soil, data = rhizo_chao_Barley_NPK_UF)

###ROOT

idx <- index_root$Soil %in% c("NPK","PK")
root_chao_Barley_NPK_PK <- index_root[idx, ]

wilcox.test(value~Soil, data = root_chao_Barley_NPK_PK)

idx <- index_root$Soil %in% c("PK","UF")
root_chao_Barley_PK_UF <- index_root[idx, ]

wilcox.test(value~Soil, data = root_chao_Barley_PK_UF)

idx <- index_root$Soil %in% c("NPK","UF")
root_chao_Barley_NPK_UF <- index_root[idx, ]

wilcox.test(value~Soil, data = root_chao_Barley_NPK_UF)


###NODULES

idx <- index_nodules$Soil %in% c("NPK","PK")
nodules_chao_Barley_NPK_PK <- index_nodules[idx, ]

wilcox.test(value~Soil, data = nodules_chao_Barley_NPK_PK)

idx <- index_nodules$Soil %in% c("PK","UF")
nodules_chao_Barley_PK_UF <- index_nodules[idx, ]

wilcox.test(value~Soil, data = nodules_chao_Barley_PK_UF)

idx <- index_nodules$Soil %in% c("NPK","UF")
nodules_chao_Barley_NPK_UF <- index_nodules[idx, ]

wilcox.test(value~Soil, data = nodules_chao_Barley_NPK_UF)

# statistical analysis ANOVA
library(car)
library(multcompView)

index_rhizo$Soil <- as.character(index_rhizo$Soil)
index_rhizo$Compartment <- as.character(index_rhizo$Compartment)
index_rhizo$Genotype <- as.character(index_rhizo$Genotype)

ano <- aov(value ~ Soil, data=index_rhizo)

anova(ano) ### shows significant difference

### Multiple pairwise-comparsions use pair-wise t test

pairwise <- TukeyHSD(ano)


### Generate lables for significance


generate_label_df <- function(pairwise, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the plot :
  ###Tukey.labels$Type=rownames(Tukey.labels)
  ###Tukey.labels=Tukey.labels[order(Tukey.labels$Type) , ]
  return(Tukey.labels)
}


#delete rows that comtain NA values in data

pairwise[["Soil"]] <- na.omit(pairwise[["Soil"]])

# Generate the significant labels for each of my sample
## This label can be used for adding labels on plots
LABELS=generate_label_df(pairwise , "Soil")

## output Tukey results and label results

write.csv(LABELS, file =  "Askov_Lotus_rfd_chao1_ANOVA_TukeyHSD_rhizo.csv")

index_root$Soil <- as.character(index_root$Soil)
index_root$Compartment <- as.character(index_root$Compartment)
index_root$Genotype <- as.character(index_root$Genotype)

ano <- aov(value ~ Soil, data=index_root)

anova(ano) ### shows significant difference

### Multiple pairwise-comparsions use pair-wise t test

pairwise <- TukeyHSD(ano)


### Generate lables for significance


generate_label_df <- function(pairwise, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the plot :
  ###Tukey.labels$Type=rownames(Tukey.labels)
  ###Tukey.labels=Tukey.labels[order(Tukey.labels$Type) , ]
  return(Tukey.labels)
}


#delete rows that comtain NA values in data

pairwise[["Soil"]] <- na.omit(pairwise[["Soil"]])

# Generate the significant labels for each of my sample
## This label can be used for adding labels on plots
LABELS=generate_label_df(pairwise , "Soil")

## output Tukey results and label results

write.csv(LABELS, file =  "Askov_Lotus_rfd_chao1_ANOVA_TukeyHSD_root.csv")

index_nodules$Soil <- as.character(index_nodules$Soil)
index_nodules$Compartment <- as.character(index_nodules$Compartment)
index_nodules$Genotype <- as.character(index_nodules$Genotype)

ano <- aov(value ~ Soil, data=index_nodules)

anova(ano) ### shows significant difference

### Multiple pairwise-comparsions use pair-wise t test

pairwise <- TukeyHSD(ano)


### Generate lables for significance


generate_label_df <- function(pairwise, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the plot :
  ###Tukey.labels$Type=rownames(Tukey.labels)
  ###Tukey.labels=Tukey.labels[order(Tukey.labels$Type) , ]
  return(Tukey.labels)
}


#delete rows that comtain NA values in data

pairwise[["Soil"]] <- na.omit(pairwise[["Soil"]])

# Generate the significant labels for each of my sample
## This label can be used for adding labels on plots
LABELS=generate_label_df(pairwise , "Soil")

## output Tukey results and label results

write.csv(LABELS, file =  "Askov_Lotus_rfd_chao1_ANOVA_TukeyHSD_nodules.csv")