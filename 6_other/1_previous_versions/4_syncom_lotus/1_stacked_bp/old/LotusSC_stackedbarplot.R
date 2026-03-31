options(warn=-1)

# cleanup
rm(list=ls())

# load packages
library(dplyr)

# directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# files
design.file <- paste("LotusCSSP_LjSC_metadata.txt", sep="")
taxonomy.file <- paste("LjSC_taxonomy.txt", sep="")
otu_table.file <- paste("LotusCSSP_LjSC_ASVtable.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names=1, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T)

# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

idx <- rownames(otu_table) %in% taxonomy[,1]
otu_table <- otu_table[idx, ]

idx <- match(rownames(otu_table), taxonomy[,1])
taxonomy <- taxonomy[idx, ]

# otu_table normalization
otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))
df <- as.data.frame(otu_table_norm)

# add taxonomy info to "df"
tax_order <- taxonomy[ ,5]
df$order <- tax_order

ASVid <- row.names(df)
df <- cbind(ASVid, df)
row.names(df) <- NULL

### reshape data: the aim of this step is to transfer wide data to long data(check online the difference of wide and long dataset in r)
df.long <- reshape(df, varying = 2:60, direction = "long", idvar = 'ASVid', timevar="sampleID", v.names=c("RA"), 
times=c(colnames(df[,2:60])), sep="")

row.names(df.long) <- NULL

df.long <- as.data.frame(df.long)

# the following is to make stackbar plot, the top 10 order are marked with 10 different color, the rest orders are assigned to the same color

###  make a dataframe that sum RA by tax and sampleID
df.long_order <- df.long %>% select (2:4) ### select the column 2,3,4 from df.long dataframe

df.long_order <- df.long_order %>%
  group_by_(.dots = c("order", "sampleID")) %>%   ### this is to summarize RA by both tax and sampleID
  summarise_all(funs(sum))

### Add genotype information for df.long (because the OTU_table you send me is specific for the endosphere compartment and in only one soil condition, so here we only add genotype info)
df2 <- data.frame(sampleID = design$SampleID, Genotype=design$Genotype, Compartment = design$Compartment)

df.long_order <- full_join(df.long_order, df2, by="sampleID")


idx <- df.long_order$Genotype%in% c("WT","symrk","ccamk","nsp1","nsp2")
df.long_order_nof6h1 <- df.long_order[idx,]


idx <- df.long_order_nof6h1$Compartment%in% c("rhizosphere")
df.long_order_rhizo <- df.long_order_nof6h1[idx,]

idx <- df.long_order_nof6h1$Compartment%in% c("root")
df.long_order_endo <- df.long_order_nof6h1[idx,]

idx <- df.long_order_nof6h1$Compartment%in% c("nodules")
df.long_order_nod <- df.long_order_nof6h1[idx,]


## select the top 10 taxa by relative abundance for the stack bar plot

### calculate the mean RA for each taxa across samples

##Lj_CSSP
mean_RA_order_nof6h1 <- df.long_order_nof6h1 %>%
  group_by(order) %>%
  dplyr::summarize(Mean = mean(RA, na.rm = TRUE))

mean_RA_order_nof6h1 <- as.data.frame(mean_RA_order_nof6h1)

mean_RA_order_nof6h1 <- mean_RA_order_nof6h1[order(-mean_RA_order_nof6h1$Mean),]

mean_RA_order_nof6h1 %>% top_n(20)


### replace the name of tax that are not included in the top10 
library(stringr)
unique(mean_RA_order_nof6h1$order)

# stacked bar plot for the top10 taxa
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggh4x)

### The way for facel order on the plot
df.long_order_nof6h1$Genotype <- factor(df.long_order_nof6h1$Genotype, levels = c("WT","symrk","ccamk","nsp1","nsp2"))
df.long_order_nof6h1$Compartment <- factor(df.long_order_nof6h1$Compartment, levels = c("rhizosphere","root","nodules"))

### visualization by ggplot2
main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))


levels(df.long_order_nof6h1$Genotype)<- c("WT","symrk","ccamk","nsp1","nsp2") ### This is to make the letter italic for the genotype name

colors <- data.frame(group=c("Actinomycetales","Bacillales","Burkholderiales","Caulobacterales","Flavobacteriales","Pseudomonadales","Rhizobiales","Sphingomonadales","Xanthomonadales"), 
                      colors=c("#2A0134","#771155", "#645394", "#AA4488","#ffeeef", "#88CCAA", "lightyellow", "#fdbb6b", "#ffc0cb")) ### color used here is from "Paired"
# 
colors <- colors[colors$group %in%df.long_order_nof6h1$order, ]

# colors <- data.frame(group=c("Burkholderiales","Caulobacterales","Chloroflexales","Flavobacteriales","Frankiales","Gaiellales","Gemmatimonadales","Micrococcales","Micromonosporales","Propionibacteriales","Pseudomonadales","Pseudonocardiales","Rhizobiales","Sphingomonadales","Streptomycetales","Unknown","Xanthomonadales","Other"), 
#                      colors=c("#645394","#AA4488","#CC99BB","#ffeeef","#114477","#4477AA","#77AADD","#44AAAA","#77CCCC","#117744","#88CCAA","#CDEBC5","lightyellow","#fdbb6b","#ffd7b5","darkgrey","#ffc0cb","lightgrey")) ### color used here is from "Paired"
# 
# colors <- colors[colors$group %in%df.long_order_nof6h1$order, ]

p1 <- ggplot(df.long_order_nof6h1, aes(x=sampleID, y = RA, fill = order)) +
  geom_bar(stat = "identity", width = .5) +
  facet_nested(~Compartment + Genotype,scales ="free_x", labeller = label_parsed)+
  scale_fill_manual(values=colors$colors)+
  main_theme+
  ylab("Relative abundance")+
  theme(legend.position = "bottom")+ guides(fill=guide_legend())+ 
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),  # Adjust the size here
        strip.text.x = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20))
p1

ggsave(paste(figures.dir, "LotusSC_stacked_barplot.png", sep=""), p1, width=14, height=8)
ggsave(paste(figures.dir, "LotusSC_stacked_barplot.pdf", sep=""), p1, width=14, height=8)


######BOXPLOTS
##both compartments
idx <- df.long_order_nof6h1$Compartment%in% c("rhizosphere","root")
df.long_order_nonod <- df.long_order_nof6h1[idx,]

colors2 <- data.frame(group=c("WT","symrk","ccamk","nsp1","nsp2"), 
                     color=c("#33a02c","#ff7f00","#1f78b4","#e31a1c", "#ffd700"))

dodge <- position_dodge (width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))
p8 <- ggplot(df.long_order_nonod, aes(x=order, y=RA, fill=Genotype)) +
  geom_boxplot( width=0.5,position = dodge, outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.3), size=1, alpha=0.5)+
  scale_fill_manual(values=as.character(colors2$color)) +
  labs(x="", y="Accumulative Relative Abundance") +
  scale_y_continuous(limits = c(0, 0.6))+
  facet_wrap(~Compartment,  ncol=1)+
  main_theme +
  theme(legend.position= "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 70, hjust=1))
p8

ggsave(paste(figures.dir, "LotusSC_boxplots.png", sep=""), p8, width=16, height=12)
ggsave(paste(figures.dir, "LotusSC_boxplots.pdf", sep=""), p8, width=16, height=12)

##rhizo only
idx <- df.long_order_nof6h1$Compartment%in% c("rhizosphere")
df.long_order_rhizo <- df.long_order_nof6h1[idx,]

colors2 <- data.frame(group=c("WT","symrk","ccamk","nsp1","nsp2"), 
                      color=c("#33a02c","#ff7f00","#1f78b4","#e31a1c", "#ffd700"))

dodge <- position_dodge (width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))
p8 <- ggplot(df.long_order_rhizo, aes(x=order, y=RA, fill=Genotype)) +
  geom_boxplot( width=1,position = dodge, outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.3), size=1, alpha=0.5)+
  scale_fill_manual(values=as.character(colors2$color)) +
  labs(x="", y="Accumulative Relative Abundance") +
  scale_y_continuous(limits = c(0, 0.5))+
  #facet_wrap(~Compartment,  ncol=1)
main_theme +
  theme(legend.position= "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 70, hjust=1))
p8

ggsave(paste(figures.dir, "LotusSC_boxplots_rhizo.png", sep=""), p8, width=20, height=10)
ggsave(paste(figures.dir, "LotusSC_boxplots_rhizo.pdf", sep=""), p8, width=20, height=10)

##root only
idx <- df.long_order_nof6h1$Compartment%in% c("root")
df.long_order_endo <- df.long_order_nof6h1[idx,]

colors2 <- data.frame(group=c("WT","symrk","ccamk","nsp1","nsp2"), 
                      color=c("#33a02c","#ff7f00","#1f78b4","#e31a1c", "#ffd700"))

dodge <- position_dodge (width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))
p8 <- ggplot(df.long_order_endo, aes(x=order, y=RA, fill=Genotype)) +
  geom_boxplot( width=1,position = dodge, outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.3), size=1, alpha=0.5)+
  scale_fill_manual(values=as.character(colors2$color)) +
  labs(x="", y="Accumulative Relative Abundance") +
  scale_y_continuous(limits = c(0, 0.5))+
  #facet_wrap(~Compartment,  ncol=1)
  main_theme +
  theme(legend.position= "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 70, hjust=1))
p8

ggsave(paste(figures.dir, "LotusSC_boxplots_endo.png", sep=""), p8, width=20, height=10)
ggsave(paste(figures.dir, "LotusSC_boxplots_endo.pdf", sep=""), p8, width=20, height=10)


#########significant ones only - BOXPLOT
idx <- df.long_order_nof6h1$Compartment%in% c("rhizosphere","root")
df.long_order_nonod <- df.long_order_nof6h1[idx,]

idx <- df.long_order_nonod$order%in% c("Burkholderiales","Flavobacteriales","Pseudomonadales","Rhizobiales")
df.long_order_nonod_significant <- df.long_order_nonod[idx,]

colors2 <- data.frame(group=c("WT","symrk","ccamk","nsp1","nsp2"), 
                      color=c("#33a02c","#ff7f00","#1f78b4","#e31a1c", "#ffd700"))

dodge <- position_dodge (width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=40),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

p8 <- ggplot(df.long_order_nonod_significant, aes(x=order, y=RA, fill=Genotype)) +
  geom_boxplot( width=0.5,position = dodge, outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.3), size=1, alpha=0.5)+
  scale_fill_manual(values=as.character(colors2$color)) +
  labs(x="", y="Accumulative Relative Abundance") +
  scale_y_continuous(limits = c(0, 0.6))+
  #facet_nested(~Compartment,scales ="free_x", labeller = label_parsed)+
  facet_wrap(~Compartment, ncol=1)+
main_theme +
  theme(legend.position= "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
p8

ggsave(paste(figures.dir, "LotusSC_boxplots_significant.png", sep=""), p8, width=16, height=12)
ggsave(paste(figures.dir, "LotusSC_boxplots_significant.pdf", sep=""), p8, width=16, height=12)

#########significant ones only - BARPLOT
main_theme <- theme(panel.background=element_blank(),
                    # panel.grid=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=15),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

p10 <- ggplot(df.long_order_nonod_significant, aes(x = order, y = RA, fill = Genotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_jitter(aes(group = Genotype), position = position_jitterdodge(jitter.width = 0.3), size = 1, alpha = 0.5) +
  scale_fill_manual(values = as.character(colors2$color)) +
  labs(x = "", y = "Cumulative Relative Abundance") +
  facet_nested(~Compartment,scales ="free_x", labeller = label_parsed)+
  # facet_wrap(~Compartment, ncol=1)+
  scale_y_continuous(limits = c(0, 0.6)) +
  main_theme +
  theme(legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 30, hjust = 1))

p10

ggsave(paste(figures.dir, "LotusSC_barplot_significant.png", sep=""), p10, width=10, height=6)
ggsave(paste(figures.dir, "LotusSC_barplot_significant.pdf", sep=""), p10, width=10, height=6)

##########SIGNIFICANCE TESTING##################
unique(df.long_order_nof6h1$order)

##################Rhizo
idx <- df.long_order_rhizo$order %in% c("Actinomycetales")
df.long_order_rhizo_Actinomycetales <- df.long_order_rhizo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_rhizo_Actinomycetales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_rhizo_Actinomycetales_ANOVA.csv")

#####
idx <- df.long_order_rhizo$order %in% c("Bacillales")
df.long_order_rhizo_Bacillales <- df.long_order_rhizo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_rhizo_Bacillales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_rhizo_Bacillales_ANOVA.csv")

#####
idx <- df.long_order_rhizo$order %in% c("Burkholderiales")
df.long_order_rhizo_Burkholderiales <- df.long_order_rhizo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_rhizo_Burkholderiales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_rhizo_Burkholderiales_ANOVA.csv")

#####
idx <- df.long_order_rhizo$order %in% c("Caulobacterales")
df.long_order_rhizo_Caulobacterales <- df.long_order_rhizo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_rhizo_Caulobacterales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_rhizo_Caulobacterales_ANOVA.csv")

#####
idx <- df.long_order_rhizo$order %in% c("Flavobacteriales")
df.long_order_rhizo_Flavobacteriales <- df.long_order_rhizo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_rhizo_Flavobacteriales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_rhizo_Flavobacteriales_ANOVA.csv")

#####
idx <- df.long_order_rhizo$order %in% c("Pseudomonadales")
df.long_order_rhizo_Pseudomonadales <- df.long_order_rhizo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_rhizo_Pseudomonadales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_rhizo_Pseudomonadales_ANOVA.csv")

#####
idx <- df.long_order_rhizo$order %in% c("Rhizobiales")
df.long_order_rhizo_Rhizobiales <- df.long_order_rhizo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_rhizo_Rhizobiales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_rhizo_Rhizobiales_ANOVA.csv")

#####
idx <- df.long_order_rhizo$order %in% c("Sphingomonadales")
df.long_order_rhizo_Sphingomonadales <- df.long_order_rhizo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_rhizo_Sphingomonadales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_rhizo_Sphingomonadales_ANOVA.csv")

#####
idx <- df.long_order_rhizo$order %in% c("Xanthomonadales")
df.long_order_rhizo_Xanthomonadales <- df.long_order_rhizo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_rhizo_Xanthomonadales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_rhizo_Xanthomonadales_ANOVA.csv")

##################endo
idx <- df.long_order_endo$order %in% c("Actinomycetales")
df.long_order_endo_Actinomycetales <- df.long_order_endo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_endo_Actinomycetales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_endo_Actinomycetales_ANOVA.csv")

#####
idx <- df.long_order_endo$order %in% c("Bacillales")
df.long_order_endo_Bacillales <- df.long_order_endo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_endo_Bacillales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_endo_Bacillales_ANOVA.csv")

#####
idx <- df.long_order_endo$order %in% c("Burkholderiales")
df.long_order_endo_Burkholderiales <- df.long_order_endo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_endo_Burkholderiales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_endo_Burkholderiales_ANOVA.csv")

#####
idx <- df.long_order_endo$order %in% c("Caulobacterales")
df.long_order_endo_Caulobacterales <- df.long_order_endo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_endo_Caulobacterales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_endo_Caulobacterales_ANOVA.csv")

#####
idx <- df.long_order_endo$order %in% c("Flavobacteriales")
df.long_order_endo_Flavobacteriales <- df.long_order_endo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_endo_Flavobacteriales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_endo_Flavobacteriales_ANOVA.csv")

#####
idx <- df.long_order_endo$order %in% c("Pseudomonadales")
df.long_order_endo_Pseudomonadales <- df.long_order_endo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_endo_Pseudomonadales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_endo_Pseudomonadales_ANOVA.csv")

#####
idx <- df.long_order_endo$order %in% c("Rhizobiales")
df.long_order_endo_Rhizobiales <- df.long_order_endo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_endo_Rhizobiales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_endo_Rhizobiales_ANOVA.csv")

#####
idx <- df.long_order_endo$order %in% c("Sphingomonadales")
df.long_order_endo_Sphingomonadales <- df.long_order_endo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_endo_Sphingomonadales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_endo_Sphingomonadales_ANOVA.csv")

#####
idx <- df.long_order_endo$order %in% c("Xanthomonadales")
df.long_order_endo_Xanthomonadales <- df.long_order_endo[idx,]
ano <- aov(RA ~ Genotype, data=df.long_order_endo_Xanthomonadales)
anova(ano)
pairwise <- TukeyHSD(ano)
library(multcompView)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")
write.csv(LABELS, file = "LotusSC_endo_Xanthomonadales_ANOVA.csv")