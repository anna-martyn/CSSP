options(warn=-1)

# cleanup

rm(list=ls())

# directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# files

weight.file <- paste("BarelyCSSP_cerealSC_diaz_shootfw.txt", sep = "")

# load data

weight <- read.table(weight.file, header=T, sep="\t")


idx <- weight$Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2")

weight_noOX <- weight[idx,]

idx <- weight_noOX$Inoculum %in% c("SC_only")

weight_noOX_SConly<- weight_noOX[idx,]


#design <- design[idx, ]

# set colour for the genotype

colors <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
                     color=c("#A9C289","#FEDA8B","#FDB366","#C0E4EF","#6EA6CD"))

shapes <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
                     shape=c(19,19,19,19,19))

# set order for plants

weight_noOX_SConly$Genotype <- factor(weight_noOX_SConly$Genotype, levels = colors$group)
weight_noOX_SConly$Genotype <- factor(weight_noOX_SConly$Genotype, levels = shapes$group)

# load plotting functions

###install.packages("ggplot2")

library("ggplot2")
library("scales")
library("grid")
library(RColorBrewer)

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

p2 <- ggplot(weight_noOX_SConly, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  geom_boxplot(aes(shape=Genotype), width=0.3,position = dodge,outlier.color = NA)+
  geom_jitter(aes(shape=Genotype), position=position_jitterdodge(jitter.width = 0.3), size=3, alpha=0.3)+
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


p2

ggsave(paste("BarleySConly_shootfw.png", sep=""), p2, width=4, height=6)
ggsave(paste("BarleySConly_shootfw.pdf", sep=""), p2, width=4, height=6)



library(car)

ano <- aov(Fresh_weight ~ Genotype, data=weight_noOX_SConly)

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

write.csv(LABELS, file =  "Fresh_weights_TukeyHSD_BarleySConly.csv")



