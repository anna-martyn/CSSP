# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data.
design <- read.table("BarleyCSSP_Askov_reseq_metadata.txt", header=T, sep="\t")
taxonomy <- read.table("Barley_Askov_Rep_10_4_taxonomy.txt", sep="\t", header=T, fill=T)
asv_table <- read.table("BarleyCSSP_Askov_reseq_ASVtable_10_4.txt", sep="\t", header=T, row.names =1, check.names=F)

# Load required packages.
library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

# Clean-up taxonomy file layout.
taxonomy <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right") %>%
  mutate(across(Kingdom:Species, ~sub("^.{3}", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Remove the three Lotus soil samples from the dataset, and only keep design samples which are present in ASV table.
design <- design %>%
  filter(!(Genotype == "Soil" & grepl("Lj", Description))) %>%
  filter(Sample_ID %in% colnames(asv_table))

# Subset and reorder ASV table to match filtered design.
asv_table <- asv_table %>%
  select(all_of(design$Sample_ID)) %>%
  filter(rownames(.) %in% taxonomy$ASVid)

# Subset dataset to only keep soil samples.
design_bulk <- design %>% filter(Genotype == "Soil")
asv_table_bulk <- asv_table[, design_bulk$Sample_ID]

# Convert ASV reads to relative abundances.
asv_table_norm <- apply(asv_table_bulk, 2, function(x) x / sum(x))

# Calculate beta-diversity using Bray-Curtis distances.
bray_curtis <- vegdist(t(asv_table_norm), method="bray")

# Make PCoA plot.

## Define colours and shapes for graph.
# colors <- data.frame(group=c("NPK","PK","UF"), 
#                      colors=c("#341C02","#A06A37","#D2B48C"))

# colors <- data.frame(group=c("NPK","PK","UF"),
#                      colors=c("#73675A","#A89F90","#D6C9B8"))

# colors <- data.frame(group=c("NPK","PK","UF"),
#                      colors=c("#9CAF88","#D08970","#7A9E9F"))

# colors <- data.frame(group=c("NPK","PK","UF"),
#                      colors=c("#A8D17F","#E07B5F","#5FB0B7"))

colors <- data.frame(group=c("NPK","PK","UF"),
                     colors=c("#6F944F","#B2563C","#3C7D82"))

# colors <- data.frame(group=c("NPK","PK","UF"),
#                      colors=c("#A89C94","#B2A083","#9D7F6D"))

# colors <- data.frame(group=c("NPK","PK","UF"), 
#                      colors=c("#83695C","#B89A86","#E2CFC1"))

# colors <- data.frame(group=c("NPK","PK","UF"), 
#                      colors=c("#6D3B00","#A06A37","#D2B48C"))

shapes <- data.frame(group=c("NPK","PK", "UF"),
                     shape=c(19, 19, 19))

colors <- colors[colors$group %in% design_bulk$Soil, ]
shapes <- shapes[shapes$group %in% design_bulk$Soil,]


# Perform principal coordinate analysis of Bray-Curtis distances.
k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design_bulk[match(rownames(points), design_bulk$Sample_ID), ])

colors <- colors[colors$group %in% points$Soil, ]
points$Soil <- factor(points$Soil, levels=colors$group)

shapes <- shapes[shapes$group %in% points$Soil, ]
points$Soil <- factor(points$Soil, levels=shapes$group)

# calculate centroids per Soil group
centroids <- aggregate(cbind(points$x, points$y) ~ Soil, data=points, FUN=mean)

# join centroids back to the points so each sample has its group centroid
segments <- merge(points, setNames(centroids, c('Soil','seg_x','seg_y')), 
                  by="Soil", sort=FALSE)

# Plot PCo 1 and 2.
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 8, color = "black"),
                    legend.text = element_text(size=8, color = "black"),
                    legend.title = element_text(size=8, color = "black",
                                                hjust = 0.5, margin = margin(b = -2)),
                    legend.key=element_blank(),
                    axis.title.y = element_text(size = 8),
                    text=element_text(size=8, color="black"),
                    legend.position="right",
                    legend.margin = margin(l = -10),
                    # legend.background=element_rect(colour="black", fill=NA),
                    legend.background=element_blank(),
                    plot.title = element_text(size=8, hjust=0.9),
                    legend.key.spacing.y = unit(-0.3, 'cm'))

# Old format with browns
# p <- ggplot(points, aes(x=x, y=y, color=Soil, shape=Soil)) +
#   geom_point(alpha=0.7, size=2) +
#   geom_segment(data=segments,
#                mapping=aes(x=x, y=y, xend=seg_x, yend=seg_y, color=Soil),
#                alpha=0.5) +
#   scale_colour_manual(values=as.character(colors$colors)) +
#   scale_shape_manual(values=shapes$shape) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
#   main_theme

p <- ggplot(points, aes(x=x, y=y, fill=Soil)) +
  geom_point(alpha=0.7, size=3, shape = 21) +
  geom_segment(data=segments,
               mapping=aes(x=x, y=y, xend=seg_x, yend=seg_y, fill=Soil),
               alpha=0.5) +
  scale_fill_manual(values=as.character(colors$colors)) +
  # scale_shape_manual(values=shapes$shape) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  main_theme

p

# Save plot.
ggsave(paste("Barley_bulk_PCoA.pdf", sep=""),
       p, width=5, height=5, units = "cm")
saveRDS(p, file = "Barley_bulk_PCoA.rds")
saveRDS(p, file = "../6_final_figure/Barley_bulk_PCoA.rds")

