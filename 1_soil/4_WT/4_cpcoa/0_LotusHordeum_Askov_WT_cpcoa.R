# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
library("ggplot2")
library("vegan")

# Load the cpcoa functions file.
source("cpcoa.func.R")

# Load the Lotus files.
Lotus_design <- read.table(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", header=T, sep="\t"
)
Lotus_asv_table <- read.table(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t", header = TRUE, row.names = 1, 
  check.names = FALSE, comment.char = ""
)

# Load Hordeum files.
Hordeum_design <- read.table(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt", header=T, sep="\t"
)
Hordeum_asv_table <- read.table(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t", header = TRUE, row.names = 1, 
  check.names = FALSE, comment.char = "", skip = 1
)

# Keep the samples that are present in the design and ASV table files for each dataset.
idx<- intersect(Lotus_design$SampleID, colnames(Lotus_asv_table))
Lotus_design <- Lotus_design[match(idx, Lotus_design$SampleID), ]
Lotus_asv_table <- Lotus_asv_table[, idx]

idx<- intersect(Hordeum_design$SampleID, colnames(Hordeum_asv_table))
Hordeum_design <- Hordeum_design[match(idx, Hordeum_design$SampleID), ]
Hordeum_asv_table <- Hordeum_asv_table[, idx]

# Subset the datasets to only keep samples with genotype WT.
Lotus_wt_samples <- Lotus_design$SampleID[Lotus_design$Genotype == "WT"]
Lotus_design <- Lotus_design[Lotus_design$SampleID %in% Lotus_wt_samples, ]
Lotus_asv_table <- Lotus_asv_table[, Lotus_wt_samples]

Hordeum_wt_samples <- Hordeum_design$SampleID[Hordeum_design$Genotype == "WT"]
Hordeum_design <- Hordeum_design[Hordeum_design$SampleID %in% Hordeum_wt_samples, ]
Hordeum_asv_table <- Hordeum_asv_table[, Hordeum_wt_samples]

# Calculate the Bray-Curtis distances.
Lotus_asv_table_norm <- apply(Lotus_asv_table, 2, function(x) x / sum(x))
bray_curtis_Lotus <- vegdist(t(Lotus_asv_table_norm), method="bray")

Hordeum_asv_table_norm <- apply(Hordeum_asv_table, 2, function(x) x / sum(x))
bray_curtis_Hordeum <- vegdist(t(Hordeum_asv_table_norm), method="bray")

# Set the colours and shapes for the plot.
colors <- data.frame(group=c("NPK","PK","UF"),
                     colors=c("#6F944F","#B2563C","#3C7D82"))

shapes <- data.frame(group=c("Rhizosphere","Root","Nodules"), shapes=c(15,16,17))

# Set the main theme for the plots.
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA,
                                                linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 8, color = "black"),
                    legend.text = element_text(size = 8, color = "black"),
                    legend.key=element_blank(),
                    axis.title.y = element_text(size = 8),
                    text=element_text(size = 8, color="black"),
                    legend.position="none",
                    legend.background=element_blank(),
                    plot.title = element_text(size = 8, hjust=0.9))

# Perform the Lotus cpcoa analysis.
sqrt_transform <- T
capscale.gen <- capscale(
  t(Lotus_asv_table_norm) ~ Compartment*Soil,
  data = Lotus_design,
  add = F, 
  sqrt.dist = sqrt_transform, 
  distance = "bray"
)

## Perform an ANOVA-like permutation analysis.
perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

## Generate variability tables and calculate the confidence intervals for the variance.
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

## Calculate the centroids per soil-compartment combination.
centroids <- aggregate(cbind(x,y) ~ Soil + Compartment, data=points, FUN=mean)
segments <- merge(points, setNames(centroids, c("Soil","Compartment","seg_x","seg_y")),
                  by=c("Soil","Compartment"), sort=FALSE)

## Plot CPCo 1 and 2.
p1 <- ggplot(points, aes(x=x, y=y, color=Soil, shape=Compartment)) +
  geom_segment(data=segments, aes(x=x, y=y, xend=seg_x, yend=seg_y, color=Soil),
               alpha=0.5, show.legend=FALSE) +
  geom_point(size=3, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("CPCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  ggtitle("Lotus", subtitle = paste(format(100 * variance, digits=3), 
                                    "% of variance; p<", format(p.val, digits=2), sep="")) +
  main_theme +
  theme(legend.position="right",
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        plot.title = element_text(face="bold", size = 8, hjust=0),
        plot.subtitle = element_text(size = 8, hjust=0)) 

p1

## Save the plot.
ggsave(paste("Lotus_Askov_WT_cpcoa.pdf", sep=""), p1, width=5, height=5, units = "cm")
saveRDS(p1, file = "Lotus_Askov_WT_cpcoa.rds")
saveRDS(p1, file = "../7_final_figures/Lotus_Askov_WT_cpcoa.rds")

# Perform the Hordeum cpcoa analysis.
sqrt_transform <- T
capscale.gen <- capscale(
  t(Hordeum_asv_table_norm) ~ Compartment*Soil, 
  data = Hordeum_design, 
  add = F, 
  sqrt.dist = sqrt_transform, 
  distance = "bray"
)

## Perform an ANOVA-like permutation analysis.
perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

## Generate variability tables and calculate the confidence intervals for the variance.
var_tbl.gen <- variability_table(capscale.gen)
eig <- capscale.gen$CCA$eig
variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

## Extract the weighted average (sample) scores.
points <- capscale.gen$CCA$wa[, 1:2]
points <- as.data.frame(points)

colnames(points) <- c("x", "y")

points <- cbind(points, Hordeum_design[match(rownames(points), Hordeum_design$SampleID), ])

points$Soil <- factor(points$Soil, levels=colors$group)
points$Compartment <- factor(points$Compartment, levels=shapes$group)

## Calculate the centroids per soil-compartment combination.
centroids <- aggregate(cbind(x,y) ~ Soil + Compartment, data=points, FUN=mean)
segments <- merge(points, setNames(centroids, c("Soil","Compartment","seg_x","seg_y")),
                  by=c("Soil","Compartment"), sort=FALSE)


## Plot CPCo 1 and 2.
p2 <- ggplot(points, aes(x=x, y=y, color=Soil, shape=Compartment)) +
  geom_segment(data=segments, aes(x=x, y=y, xend=seg_x, yend=seg_y, color=Soil),
               alpha=0.5, show.legend=FALSE) +
  geom_point(size=3, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  labs(x=paste("CPCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  ggtitle("Hordeum", subtitle = paste(format(100 * variance, digits=3), 
                                    "% of variance; p<", format(p.val, digits=2), sep="")) +
  main_theme +
  theme(legend.position="right",
        legend.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        plot.title = element_text(face="bold", size = 8, hjust=0),
        plot.subtitle = element_text(size = 8, hjust=0)) 

p2

## Save the plot.
ggsave(paste("Hordeum_Askov_WT_cpcoa.pdf", sep=""), p2, width=5, height=5)
saveRDS(p2, file = "Hordeum_Askov_WT_cpcoa.rds")
saveRDS(p2, file = "../7_final_figures/Hordeum_Askov_WT_cpcoa.rds")
