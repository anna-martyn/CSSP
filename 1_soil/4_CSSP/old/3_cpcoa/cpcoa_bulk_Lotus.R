
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de
#

options(warn=-1)

# cleanup

rm(list=ls())

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load plotting functions
library("ggplot2")
library("scales")
library("grid")
library("vegan")

# load plotting functions directory
source("cpcoa.func.R")

# files
design.file <- paste("Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt", sep="")
taxonomy.file <- paste("LotusSep_exclUFnew_10_4_silva138_taxonomy.txt", sep="")
otu_table.file <- paste("LotusSep_exclUFnew_ASVtable_10_4_nospike.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names =1, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T)

# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

idx <- rownames(otu_table) %in% taxonomy[, 1]
otu_table <- otu_table[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# subset samples by fraction
idx <- design$Compartment %in% c( "rhizosphere")
design_rhizo <- design[idx, ]
otu_table_rhizo <- otu_table[, idx]

idx <- design$Compartment %in% c( "root")
design_root <- design[idx, ]
otu_table_root <- otu_table[, idx]

# calculate Bray-Curtis distances
otu_table_norm <- apply(otu_table_rhizo, 2, function(x) x / sum(x))
bray_curtis <- vegdist(t(otu_table_norm), method="bray")

# CPCoA
# colors <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
#                      color=c("#33a02c","#ff7f00","#1f78b4","#e31a1c","#ffd700"))
colors <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
                     color=c("#A9C289","#FEDA8B","#FDB366","#C0E4EF","#6EA6CD"))
colors <- colors[colors$group %in% design_rhizo$Genotype, ]

shapes <- data.frame(group=c("NPK","PK", "UF"),
                     shape=c(18, 17, 16))
shapes <- shapes[shapes$group %in% design_rhizo$Soil, ]


sqrt_transform <- T


capscale.gen <- capscale(t(otu_table_norm) ~ Genotype*Soil + Condition(Bio_rep), data=design_rhizo, add=F, sqrt.dist=sqrt_transform, distance="bray")


# ANOVA-like permutation analysis

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

# generate variability tables and calculate confidence intervals for the variance

var_tbl.gen <- variability_table(capscale.gen)

eig <- capscale.gen$CCA$eig

variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

# extract the weighted average (sample) scores

points <- capscale.gen$CCA$wa[, 1:2]
points <- as.data.frame(points)

colnames(points) <- c("x", "y")

points <- cbind(points, design_rhizo[match(rownames(points), design_rhizo$SampleID), ])

points$Genotype <- factor(points$Genotype, levels=colors$group)
points$Soil <- factor(points$Soil, levels=shapes$group)


# plot CPCo 1 and 2
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line.x=element_line(color="black", size=1),
                    axis.line.y=element_line(color="black", size=1),
                    axis.ticks=element_line(color="black", size=1),
                    axis.text=element_text(colour="black", size=30),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans",, color="black"))

p1 <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Soil)) +
  stat_ellipse(type = "norm", level = 0.8)+
  geom_point(size=8) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape)+
  labs(x=paste("CPCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
                format(p.val, digits=2),
                sep="")) +
  main_theme +
  theme(legend.position="right", 
        plot.title = element_text(size = 30, face="bold"), 
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

p1

ggsave(paste("cpcoa_Genotype_Soil_Lotus_rhizo_NEW.png", sep=""), p1, width=10, height=8)
ggsave(paste("cpcoa_Genotype_Soil_Lotus_rhizo_NEW.pdf", sep=""), p1, width=10, height=8)

# calculate Bray-Curtis distances
otu_table_norm <- apply(otu_table_root, 2, function(x) x / sum(x))
bray_curtis <- vegdist(t(otu_table_norm), method="bray")

# CPCoA
# colors <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
#                      color=c("#33a02c","#ff7f00","#1f78b4","#e31a1c","#ffd700"))
colors <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
                     color=c("#A9C289","#FEDA8B","#FDB366","#C0E4EF","#6EA6CD"))
colors <- colors[colors$group %in% design_root$Genotype, ]

shapes <- data.frame(group=c("NPK","PK", "UF"),
                     shape=c(18, 17, 16))
shapes <- shapes[shapes$group %in% design_root$Soil, ]


sqrt_transform <- T


capscale.gen <- capscale(t(otu_table_norm) ~ Genotype*Soil + Condition(Bio_rep), data=design_root, add=F, sqrt.dist=sqrt_transform, distance="bray")


# ANOVA-like permutation analysis

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

# generate variability tables and calculate confidence intervals for the variance

var_tbl.gen <- variability_table(capscale.gen)

eig <- capscale.gen$CCA$eig

variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

# extract the weighted average (sample) scores

points <- capscale.gen$CCA$wa[, 1:2]
points <- as.data.frame(points)

colnames(points) <- c("x", "y")

points <- cbind(points, design_root[match(rownames(points), design_root$SampleID), ])

points$Genotype <- factor(points$Genotype, levels=colors$group)
points$Soil <- factor(points$Soil, levels=shapes$group)


# plot CPCo 1 and 2
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line.x=element_line(color="black", size=1),
                    axis.line.y=element_line(color="black", size=1),
                    axis.ticks=element_line(color="black", size=1),
                    axis.text=element_text(colour="black", size=30),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans",, color="black"))

p2 <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Soil)) +
  stat_ellipse(type = "norm", level = 0.8)+
  geom_point(size=8) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape)+
  labs(x=paste("CPCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
                format(p.val, digits=2),
                sep="")) +
  main_theme +
  theme(legend.position="right", 
        plot.title = element_text(size = 30, face="bold"), 
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

p2

ggsave(paste("cpcoa_Genotype_Soil_Lotus_root_NEW.png", sep=""), p2, width=10, height=8)
ggsave(paste("cpcoa_Genotype_Soil_Lotus_root_NEW.pdf", sep=""), p2, width=10, height=8)


