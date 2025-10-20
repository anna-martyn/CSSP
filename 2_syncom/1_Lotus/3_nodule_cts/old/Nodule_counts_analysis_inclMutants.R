
options(warn=-1)

# cleanup

rm(list=ls())

# directory
setwd("/Users/amartyn/Desktop/LotusSC/1_shootfw_nodcts/")
results.dir <- "/Users/amartyn/Desktop/LotusSC/1_shootfw_nodcts/"
figures.dir <- "/Users/amartyn/Desktop/LotusSC/1_shootfw_nodcts/"

# files

nod.file <- paste(results.dir, "Input_nodule_counts_inclMutants.txt", sep = "")

# load data

nod <- read.table(nod.file, header=T, sep="\t")

nod <- as.data.frame(nod)

nod$X <- NULL



idx <- nod$Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2")

nod2 <- nod[idx, ]


# set color for the soil_type

colors <- data.frame(group=c("pink", "white"), 
                     color=c("pink","lightgrey")) ### color used here is from "Paired"

shapes <- data.frame(group=c("pink","white"), shape=c(19,19))
l2 <- c("pink","white")
nod2$Nodule_type <- factor(nod2$Nodule_type, levels = l2)
shapes <- shapes[match(l2, shapes$group),]

# set order for plants

nod2$Nodule_type <- factor(nod2$Nodule_type, levels = colors$group)

nod2$Genotype <- factor(nod2$Genotype, levels = c("WT","symrk","ccamk","nsp1","nsp2"))

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
#                     axis.text=element_text(colour="black", size=20),
#                     legend.position=c(.85,.7),
#                     legend.background=element_rect(linetype=1, size=1, color = "black"),
#                     legend.key=element_blank(),
#                     text=element_text(family="sans"))

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1),
                    axis.text.y = element_text(size = 20),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))


## boxplot

dodge <- position_dodge(width = 0.75)

p1 <- ggplot(nod2, aes(x=Genotype, y=Number, fill=Nodule_type)) +
  geom_boxplot(width=0.5,position = dodge,outlier.color = NA)+
  geom_jitter(aes(shape=Genotype), position=position_jitterdodge(jitter.width = 0.3), size=3, alpha=0.3)+
  # geom_jitter(aes(), position=position_jitterdodge(jitter.width = 0.2), size=1)+
  # #facet_wrap(~Treatment, scales = "free_x", nrow = 1)+
  scale_fill_manual(values=as.character(colors$color)) +
  main_theme +
  ylab("Nodule counts/plant")+
  theme(legend.position="right", 
        plot.title = element_text(size = 20, face="bold"),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1, colour = "black"),
        legend.key.size = unit(1,"cm"))
  # theme(legend.position=c(.7,.7), 
  #       plot.title = element_text(size = 20, face="bold"), 
  #       legend.title = element_text(size = 20),
  #       strip.text.x = element_text(size = 20),
  #       legend.text = element_text(size = 20),
  #       axis.title.x = element_blank(),
  #       axis.title.y = element_text(size = 20),
  #       axis.text.x = element_text(size=20,angle=70, hjust=1))


p1

ggsave(paste(figures.dir, "LotusCSSP_LjSC_nodule_counts.png", sep=""), p1, width=6, height=6)
ggsave(paste(figures.dir, "LotusCSSP_LjSC_nodule_counts.pdf", sep=""), p1, width=6, height=6)
