options(warn=-1)

# cleanup
rm(list=ls())

# directory
setwd("/Users/amartyn/Desktop/LotusSC/1_shootfw_nodcts/")
results.dir <- "/Users/amartyn/Desktop/LotusSC/1_shootfw_nodcts/"
figures.dir <- "/Users/amartyn/Desktop/LotusSC/1_shootfw_nodcts/"

# files

weight.file <- paste(results.dir, "Input_shoot_freshweights_harvest_magentas.txt", sep = "")

# load data

weight <- read.table(weight.file, header=T, sep="\t")
weight <- read.table("Input_shoot_freshweights_harvest_magentas.txt", header=T, sep="\t")



weight <- as.data.frame(weight)
weight$X <- NULL

weight$Fresh_weight <- gsub("\\,",".",weight$Fresh_weight)

weight$Fresh_weight <- as.numeric(weight$Fresh_weight)


idx <- weight$Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2")

weight_nof6h1<- weight[idx,]

idx <- weight_nof6h1$Treatment %in% c("Lj_SC")

weight_nof6h1_SC<- weight_nof6h1[idx,]


#design <- design[idx, ]

# set colour for the genotype

colors <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
                     color=c("#33a02c","#ff7f00","#1f78b4","#e31a1c", "#ffd700"))

shapes <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
                     shape=c(19,19,19,19,19))

# set order for plants

weight_nof6h1_SC$Genotype <- factor(weight_nof6h1_SC$Genotype, levels = colors$group)
weight_nof6h1_SC$Genotype <- factor(weight_nof6h1_SC$Genotype, levels = shapes$group)

# load plotting functions

###install.packages("ggplot2")

library("ggplot2")
library("scales")
library("grid")
library(RColorBrewer)

# main_theme <- theme(panel.background=element_blank(),
#                     panel.grid=element_blank(),
#                     axis.line.x=element_line(color="black"),
#                     axis.line.y=element_line(color="black"),
#                     axis.ticks=element_line(color="black"),
#                     axis.text=element_text(colour="black", size=30),
#                     legend.position="none",
#                     legend.background=element_blank(),
#                     legend.key=element_blank(),
#                     text=element_text(family="sans"))

dodge <- position_dodge(width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text.x = element_text(size = 25, angle = 45, vjust = 1, hjust=1),
                    axis.text.y = element_text(size = 25),
                    #legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

## violin plot

dodge <- position_dodge(width = 0.4)


p2 <- ggplot(weight_nof6h1_SC, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  #geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  # geom_boxplot(width=0.4,position = dodge,outlier.color = NA)+
  geom_boxplot(aes(shape=Genotype), width=0.3,position = dodge,outlier.color = NA)+
  geom_jitter(aes(shape=Genotype), position=position_jitterdodge(jitter.width = 0.3), size=3, alpha=0.3)+
  # geom_jitter(aes(), position=position_jitterdodge(jitter.width = 0.3), size=2, alpha=0.3)+
  scale_fill_manual(values=as.character(colors$color)) +
  main_theme +
  ylab("Shoot fresh weight/plant(g)")+
  #scale_y_continuous(breaks=c(0,0.02,0.04,0.06,0.08,1.0))+
  scale_y_continuous()+
  theme(legend.position="none", 
        plot.title = element_text(size = 25, face="bold"),
        legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25),
        legend.text = element_text(size = 25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 25, colour = "black"),
        axis.text.x = element_text(size = 25, angle = 45, vjust = 1, hjust=1, colour = "black"),
        legend.key.size = unit(1,"cm"))
  # theme(legend.position="none", 
  #       plot.title = element_text(size = 30, face="bold"), 
  #       #legend.title = element_text(size = 30),
  #       strip.text.x = element_text(size = 30),
  #       legend.text = element_text(size = 30),
  #       axis.title.x = element_blank(),
  #       axis.title.y = element_text(size = 30),
  #       axis.text.x = element_text(size=30,angle = 70, hjust=1))

p2

ggsave(paste(figures.dir, "LotusSC_shootfw_nof6h1.png", sep=""), p2, width=4, height=6)
ggsave(paste(figures.dir, "LotusSC_shootfw_nof6h1.pdf", sep=""), p2, width=4, height=6)


library(car)

ano <- aov(Fresh_weight ~ Genotype, data=weight_nof6h1_SC)

anova(ano) ### shows significant difference

### Multiple pairwise-comparsions use pair-wise t test

pairwise <- TukeyHSD(ano)


### Generate lables for significance

library(multcompView)

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

pairwise[["Genotype"]] <- na.omit(pairwise[["Genotype"]])

# Generate the significant labels for each of my sample
## This label can be used for adding labels on plots
LABELS=generate_label_df(pairwise , "Genotype")

## output Tukey results and label results

write.csv(LABELS, file =  "Fresh_weights_TukeyHSD_nof6h1_SC.csv")



