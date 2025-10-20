options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Install required packages.
library("ggplot2")
library("scales")
library("grid")
library(RColorBrewer)
library(car)
library(multcompView)

# Load shoot fresh weight file.
weight.file <- paste("Input_shoot_freshweights_harvest_magentas.txt", sep = "")
weight <- read.table(weight.file, header=T, sep="\t")
weight <- read.table("Input_shoot_freshweights_harvest_magentas.txt", header=T, sep="\t")

# Modify shoot fresh weight file.
weight <- as.data.frame(weight)
weight$X <- NULL
weight$Fresh_weight <- gsub("\\,",".",weight$Fresh_weight)
weight$Fresh_weight <- as.numeric(weight$Fresh_weight)

# Select genotypes to take into account.
idx <- weight$Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2")
weight_nof6h1<- weight[idx,]

# Select only for inoculated condition.
# idx <- weight_nof6h1$Treatment %in% c("Lj_SC")
# weight_nof6h1_SC<- weight_nof6h1[idx,]

# Set colours and shapes for the genotypes
colors <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
                     color=c("#A9C289","#FEDA8B","#FDB366","#C0E4EF","#6EA6CD"))

# shapes <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
#                      shape=c(19,19,19,19,19))

shapes <- data.frame(group=c("no_SC", "Lj_SC"), 
                     shape=c(19,19))

# Set genotype orders for the plot.
weight_nof6h1$Genotype <- factor(weight_nof6h1$Genotype, levels = colors$group)
# weight_nof6h1$Genotype <- factor(weight_nof6h1$Genotype, levels = shapes$group)
weight_nof6h1$Treatment <- factor(weight_nof6h1$Treatment, levels = shapes$group)

# weight_nof6h1_SC$Genotype <- factor(weight_nof6h1_SC$Genotype, levels = colors$group)
# weight_nof6h1_SC$Genotype <- factor(weight_nof6h1_SC$Genotype, levels = shapes$group)

# Make graph.
dodge <- position_dodge(width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1),
                    axis.text.y = element_text(size = 20),
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

p2 <- ggplot(weight_nof6h1, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  geom_boxplot(aes(shape=Treatment), width=0.5,position = dodge,outlier.color = NA)+
  geom_jitter(position=position_jitterdodge(jitter.width = 0.3), size=2, alpha=0.3)+
  scale_fill_manual(values=as.character(colors$color)) +
  facet_wrap(~Treatment, scales = "free_x", nrow=1)+
  main_theme +
  ylab("Shoot fresh weight/plant [g]")+
  scale_y_continuous()+
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

ggsave(paste("LotusSC_shootfw_nof6h1_bothtreatments.png", sep=""), p2, width=4, height=6)
ggsave(paste("LotusSC_shootfw_nof6h1_bothtreatments.pdf", sep=""), p2, width=4, height=6)


# Statistics.

## Subset dataframe by treatment (uninoc WT vs. mutants, LjSC WT vs. mutants).
idx <- weight_nof6h1$Treatment %in% c("no_SC")
weight_nof6h1_uninoc<- weight_nof6h1[idx,]

idx <- weight_nof6h1$Treatment %in% c("Lj_SC")
weight_nof6h1_SC<- weight_nof6h1[idx,]

### uninoculated
ano <- aov(Fresh_weight ~ Genotype, data=weight_nof6h1_uninoc)
anova(ano)
pairwise <- TukeyHSD(ano)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}

pairwise[["Genotype"]] <- na.omit(pairwise[["Genotype"]])
LABELS=generate_label_df(pairwise , "Genotype")

write.csv(LABELS, file =  "Fresh_weights_TukeyHSD_nof6h1_uninoc.csv")

### inoculated with LjSC
ano <- aov(Fresh_weight ~ Genotype, data=weight_nof6h1_SC)
anova(ano)
pairwise <- TukeyHSD(ano)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}

pairwise[["Genotype"]] <- na.omit(pairwise[["Genotype"]])
LABELS=generate_label_df(pairwise , "Genotype")

write.csv(LABELS, file =  "Fresh_weights_TukeyHSD_nof6h1_SC.csv")

## Subset dataframe by genotype (e.g. WT uninoc. vs WT LjSC).
### WT
idx <- weight_nof6h1$Genotype %in% c("WT")
weight_nof6h1_WT<- weight_nof6h1[idx,]

pairwise.wilcox.test(weight_nof6h1_WT$Fresh_weight, weight_nof6h1_WT$Treatment, p.adjust.method = "BH")

### symrk
idx <- weight_nof6h1$Genotype %in% c("symrk")
weight_nof6h1_symrk<- weight_nof6h1[idx,]

pairwise.wilcox.test(weight_nof6h1_symrk$Fresh_weight, weight_nof6h1_symrk$Treatment, p.adjust.method = "BH")

### ccamk
idx <- weight_nof6h1$Genotype %in% c("ccamk")
weight_nof6h1_ccamk<- weight_nof6h1[idx,]

pairwise.wilcox.test(weight_nof6h1_ccamk$Fresh_weight, weight_nof6h1_ccamk$Treatment, p.adjust.method = "BH")

### nsp1
idx <- weight_nof6h1$Genotype %in% c("nsp1")
weight_nof6h1_nsp1<- weight_nof6h1[idx,]

pairwise.wilcox.test(weight_nof6h1_nsp1$Fresh_weight, weight_nof6h1_nsp1$Treatment, p.adjust.method = "BH")

### nsp2
idx <- weight_nof6h1$Genotype %in% c("nsp2")
weight_nof6h1_nsp2<- weight_nof6h1[idx,]

pairwise.wilcox.test(weight_nof6h1_nsp2$Fresh_weight, weight_nof6h1_nsp2$Treatment, p.adjust.method = "BH")
