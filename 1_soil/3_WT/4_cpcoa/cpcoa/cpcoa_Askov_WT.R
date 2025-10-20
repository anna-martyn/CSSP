# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load plotting functions.
library("ggplot2")
library("vegan")

# Load cpcoa functions file.
source("cpcoa.func.R")

# Load Lotus files.
Lotus_design <- read.table("./Lotus_input/Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt", header=T, sep="\t")
Lotus_asv_table <- read.table("./Lotus_input/feature-table.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, comment.char = "", skip = 1)

# Load Hordeum files.
Hordeum_design <- read.table("./Hordeum_input/BarleyCSSP_Askov_reseq_metadata.txt", header=T, sep="\t")
Hordeum_asv_table <- read.table("./Hordeum_input/BarleyCSSP_Askov_reseq_ASVtable_10_4.txt", sep="\t", header=T, row.names =1, check.names=F)

# Keep samples that are present in design and ASV table files for each dataset.
idx<- intersect(Lotus_design$SampleID, colnames(Lotus_asv_table))
Lotus_design <- Lotus_design[match(idx, Lotus_design$SampleID), ]
Lotus_asv_table <- Lotus_asv_table[, idx]

idx<- intersect(Hordeum_design$Sample_ID, colnames(Hordeum_asv_table))
Hordeum_design <- Hordeum_design[match(idx, Hordeum_design$Sample_ID), ]
Hordeum_asv_table <- Hordeum_asv_table[, idx]

# Subset datasets to only keep samples with genotype WT.
Lotus_wt_samples <- Lotus_design$SampleID[Lotus_design$Genotype == "WT"]
Lotus_design <- Lotus_design[Lotus_design$SampleID %in% Lotus_wt_samples, ]
Lotus_asv_table <- Lotus_asv_table[, Lotus_wt_samples]

Hordeum_wt_samples <- Hordeum_design$Sample_ID[Hordeum_design$Genotype == "WT"]
Hordeum_design <- Hordeum_design[Hordeum_design$Sample_ID %in% Hordeum_wt_samples, ]
Hordeum_asv_table <- Hordeum_asv_table[, Hordeum_wt_samples]

# Change compartment names in design files.
Lotus_design$Compartment[Lotus_design$Compartment == "Endosphere/Rhizoplane"] <- "Root"
Hordeum_design$Compartment[Hordeum_design$Compartment == "rhizo"] <- "Rhizosphere"
Hordeum_design$Compartment[Hordeum_design$Compartment == "endo"] <- "Root"

# Calculate Bray-Curtis distances.
Lotus_asv_table_norm <- apply(Lotus_asv_table, 2, function(x) x / sum(x))
bray_curtis_Lotus <- vegdist(t(Lotus_asv_table_norm), method="bray")

Hordeum_asv_table_norm <- apply(Hordeum_asv_table, 2, function(x) x / sum(x))
bray_curtis_Hordeum <- vegdist(t(Hordeum_asv_table_norm), method="bray")

# Set colours and shapes for plot, as well as plotting parameters.
colors <- data.frame(group=c("NPK","PK","UF"), 
                     colors=c("#341C02","#A06A37","#D2B48C"))

shapes <- data.frame(group=c("Rhizosphere","Root","Nodules"), shapes=c(15,16,17))

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 20, color = "black"),
                    legend.text = element_text(size=20, color = "black"),
                    legend.key=element_blank(),
                    axis.title.y = element_text(size = 20),
                    text=element_text(size=20, color="black"),
                    legend.position="none",
                    # legend.background=element_rect(colour="black", fill=NA),
                    legend.background=element_blank(),
                    plot.title = element_text(size=20, hjust=0.9))

# Lotus cpcoa analysis.
sqrt_transform <- T
capscale.gen <- capscale(t(Lotus_asv_table_norm) ~ Compartment*Soil + Condition(Bio_rep), data=Lotus_design, add=F, sqrt.dist=sqrt_transform, distance="bray")

## ANOVA-like permutation analysis
perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

## Generate variability tables and calculate confidence intervals for the variance.
var_tbl.gen <- variability_table(capscale.gen)
eig <- capscale.gen$CCA$eig
variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

## Extract the weighted average (sample) scores.
points <- capscale.gen$CCA$wa[, 1:2]
points <- as.data.frame(points)

colnames(points) <- c("x", "y")

points <- cbind(points, Lotus_design[match(rownames(points), Lotus_design$SampleID), ])

points$Soil <- factor(points$Soil, levels=colors$group)
points$Compartment <- factor(points$Compartment, levels=shapes$group)

## Calculate centroids per Soil x Compartment.
centroids <- aggregate(cbind(x,y) ~ Soil + Compartment, data=points, FUN=mean)
segments <- merge(points, setNames(centroids, c("Soil","Compartment","seg_x","seg_y")),
                  by=c("Soil","Compartment"), sort=FALSE)

## Plot CPCo 1 and 2.
p1 <- ggplot(points, aes(x=x, y=y, color=Soil, shape=Compartment)) +
  geom_segment(data=segments, aes(x=x, y=y, xend=seg_x, yend=seg_y, color=Soil),
               alpha=0.5, show.legend=FALSE) +
  geom_point(size=5, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("CPCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  ggtitle("Lotus", subtitle = paste(format(100 * variance, digits=3), 
                                    "% of variance; p=", format(p.val, digits=2), sep="")) +
  main_theme +
  theme(legend.position="right",
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(face="bold", size=20, hjust=0),
        plot.subtitle = element_text(size=20, hjust=0)) 

p1

ggsave(paste("cpcoa_Lotus_WT.pdf", sep=""), p1, width=5, height=5)
# ggsave(paste("cpcoa_Lotus_WT.pdf", sep=""), p1, width=10, height=8)
saveRDS(p1, file = "cpcoa_Lotus_WT.rds")

# Hordeum cpcoa analysis.
sqrt_transform <- T
capscale.gen <- capscale(t(Hordeum_asv_table_norm) ~ Compartment*Soil + Condition(Bio_rep), data=Hordeum_design, add=F, sqrt.dist=sqrt_transform, distance="bray")

## ANOVA-like permutation analysis
perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

## Generate variability tables and calculate confidence intervals for the variance.
var_tbl.gen <- variability_table(capscale.gen)
eig <- capscale.gen$CCA$eig
variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

## Extract the weighted average (sample) scores.
points <- capscale.gen$CCA$wa[, 1:2]
points <- as.data.frame(points)

colnames(points) <- c("x", "y")

points <- cbind(points, Hordeum_design[match(rownames(points), Hordeum_design$Sample_ID), ])

points$Soil <- factor(points$Soil, levels=colors$group)
points$Compartment <- factor(points$Compartment, levels=shapes$group)

## Calculate centroids per Soil x Compartment.
centroids <- aggregate(cbind(x,y) ~ Soil + Compartment, data=points, FUN=mean)
segments <- merge(points, setNames(centroids, c("Soil","Compartment","seg_x","seg_y")),
                  by=c("Soil","Compartment"), sort=FALSE)


## Plot CPCo 1 and 2.
p2 <- ggplot(points, aes(x=x, y=y, color=Soil, shape=Compartment)) +
  geom_segment(data=segments, aes(x=x, y=y, xend=seg_x, yend=seg_y, color=Soil),
               alpha=0.5, show.legend=FALSE) +
  geom_point(size=5, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("CPCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  ggtitle("Hordeum", subtitle = paste(format(100 * variance, digits=3), 
                                    "% of variance; p=", format(p.val, digits=2), sep="")) +
  main_theme +
  theme(legend.position="right",
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(face="bold", size=20, hjust=0),
        plot.subtitle = element_text(size=20, hjust=0)) 

p2

ggsave(paste("cpcoa_Barley_WT.pdf", sep=""), p2, width=5, height=5)
# ggsave(paste("cpcoa_Barley_WT.pdf", sep=""), p2, width=8, height=10)
saveRDS(p2, file = "cpcoa_Barley_WT.rds")


