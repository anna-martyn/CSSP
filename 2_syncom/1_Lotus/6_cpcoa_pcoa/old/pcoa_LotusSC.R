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

# files
design.file <- paste("LotusCSSP_LjSC_metadata.txt", sep="")
taxonomy.file <- paste("LjSC_taxonomy.txt", sep="")
otu_table.file <- paste("LotusCSSP_LjSC_ASVtable.txt", sep="")

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

# subset by Lotus samples
idx <- design$Genotype %in% c( "WT","symrk","ccamk","nsp1","nsp2")
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]

idx <- design_subset$Compartment %in% c( "rhizosphere")
design_rhizo <- design_subset[idx, ]
otu_table_rhizo <- otu_table_subset[, idx]

idx <- design_subset$Compartment %in% c( "root")
design_endo <- design_subset[idx, ]
otu_table_endo <- otu_table_subset[, idx]

#########RHIZOSPHERE
# calculate Bray-Curtis distances
otu_table_norm <- apply(otu_table_rhizo, 2, function(x) x / sum(x))
bray_curtis <- vegdist(t(otu_table_norm), method="bray")

# CPCoA
colors <- data.frame(group=c("WT","symrk","ccamk","nsp1","nsp2"), 
                     colors=c("#A9C289","#FEDA8B","#FDB366","#C0E4EF","#6EA6CD"))

colors <- colors[colors$group %in% design_rhizo$Genotype, ]

# PCoA Bray-Curtis
bray_curtis <- vegdist(t(otu_table_norm), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design_rhizo[match(rownames(points), design_rhizo$SampleID), ])

colors <- colors[colors$group %in% points$Genotype, ]
points$Genotype <- factor(points$Genotype, levels=colors$group)

# plot PCo 1 and 2
pcoa_width <- 6
pcoa_height <- 6
pcoa_size <- 6
pcoa_alpha <- 0.7

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    panel.grid=element_blank(),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 20, color = "black"),
                    legend.background = element_rect(colour = "black", fill = NA),
                    legend.text = element_text(size=20, color = "black"),
                    legend.key=element_blank(),
                    text=element_text(size=20, color="black"),
                    legend.position="right")


p <- ggplot(points, aes(x=x, y=y, color=Genotype)) +
  geom_point(alpha=pcoa_alpha, size=pcoa_size) +
  stat_ellipse(type = "norm", level = 0.8)+
  scale_colour_manual(values=as.character(colors$color)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  #ggtitle("Rhizosphere") +
  main_theme +
  theme(legend.position="right")

p

ggsave(paste("LotusSC_rhizo_PCoA_new.jpg", sep=""), p, width=pcoa_width, height=pcoa_height)
ggsave(paste("LotusSC_rhizo_PCoA_new.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

#########ROOT
# calculate Bray-Curtis distances
otu_table_norm <- apply(otu_table_endo, 2, function(x) x / sum(x))
bray_curtis <- vegdist(t(otu_table_norm), method="bray")

# CPCoA
colors <- data.frame(group=c("WT","symrk","ccamk","nsp1","nsp2"), 
                     colors=c("#A9C289","#FEDA8B","#FDB366","#C0E4EF","#6EA6CD")) 
colors <- colors[colors$group %in% design_endo$Genotype, ]

# PCoA Bray-Curtis
bray_curtis <- vegdist(t(otu_table_norm), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design_endo[match(rownames(points), design_endo$SampleID), ])

colors <- colors[colors$group %in% points$Genotype, ]
points$Genotype <- factor(points$Genotype, levels=colors$group)

# plot PCo 1 and 2
pcoa_width <- 6
pcoa_height <- 6
pcoa_size <- 6
pcoa_alpha <- 0.7

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    panel.grid=element_blank(),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 20, color = "black"),
                    legend.background = element_rect(colour = "black", fill = NA),
                    legend.text = element_text(size=20, color = "black"),
                    legend.key=element_blank(),
                    text=element_text(size=20, color="black"),
                    legend.position="right")


p2 <- ggplot(points, aes(x=x, y=y, color=Genotype)) +
  geom_point(alpha=pcoa_alpha, size=pcoa_size) +
  stat_ellipse(type = "norm", level = 0.8)+
  scale_colour_manual(values=as.character(colors$color)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  #ggtitle("Rhizosphere") +
  main_theme +
  theme(legend.position="right")

p2

ggsave(paste("LotusSC_root_PCoA_new.jpg", sep=""), p2, width=pcoa_width, height=pcoa_height)
ggsave(paste("LotusSC_root_PCoA_new.pdf", sep=""), p2, width=pcoa_width, height=pcoa_height)

