# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages.
library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggtext)
library(cowplot)

# Load data.
design <- read.table("../1_data/without_input/HordeumSC_metadata.txt", header=T, sep="\t")
asv_table <- fread("../1_data/without_input/HordeumSC_ASVtable.tsv")
setnames(asv_table, old = names(asv_table)[1], new = "ASVid")

# Filter asv table to only keep mateched ASVs.
asv_table_matched <- asv_table[grepl("_", ASVid)]
asv_table_matched <- as.data.frame(asv_table_matched)  
rownames(asv_table_matched) <- asv_table_matched[[1]]
asv_table_matched <- asv_table_matched[, -1]

## Define colours for later graph as well as main theme. Also, indicate that mutant names shall be italic in legend.
colors <- c("WT" = "#A9C289", "symrk"="#FEDA8B", "ccamk" = "#FDB366",
            "nsp1" = "#C0E4EF", "nsp2" = "#6EA6CD")

legend_labels <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA,
                                                linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 8, color = "black"),
                    legend.text = element_text(size=8, color = "black"),
                    legend.key=element_blank(),
                    axis.title.y = element_text(size = 8),
                    text=element_text(size=8, color="black"),
                    legend.position="right",
                    # legend.background=element_rect(colour="black", fill=NA),
                    legend.background=element_blank(),
                    plot.title = element_text(size=8, hjust=0.9))

compartments <- unique(design$Compartment)
segments_list <- list()
points_list <- list()
text_list <- list()
for(i in 1:length(compartments)){
  # Subset samples
  comp <- compartments[i]
  design_sub <- design %>% filter(Compartment == comp)
  asv_sub <- asv_table_matched[, colnames(asv_table_matched) %in% design_sub$SampleID]
  
  # Normalize ASVs
  asv_norm <- apply(asv_sub, 2, function(x) x / sum(x))
  
  # Bray-Curtis & PCoA
  bray <- vegdist(t(asv_norm), method="bray")
  pcoa <- cmdscale(bray, k=2, eig=TRUE)
  points <- as.data.frame(pcoa$points)
  colnames(points) <- c("x","y")
  points <- cbind(points, design_sub[match(rownames(points), design_sub$SampleID), ])
  points$Genotype <- factor(points$Genotype, levels=names(colors))
  
  var_expl <- pcoa$eig[1:2]/sum(pcoa$eig[pcoa$eig>0])
  text_dt <- data.table(
    Compartment = comp,
    text = paste0(round(var_expl[1]*100, 1), "%", "-", 
                  round(var_expl[2]*100, 1), "%")
  )
  text_list[[i]] <- text_dt
  
  # Centroids & segments
  centroids <- aggregate(cbind(x,y) ~ Genotype, data=points, FUN=mean)
  segments <- merge(points, setNames(centroids, c("Genotype","seg_x","seg_y")),
                    by = "Genotype")
  segments_list[[i]] <- segments
  points_list[[i]] <- points
}

points <- rbindlist(points_list)
segments <- rbindlist(segments_list)
text_dt <- rbindlist(text_list)

points[,Host := "Hordeum"]
segments[,Host := "Hordeum"]
text_dt[,Host := "Hordeum"]

outdir <- "./additional_outputs"

# Create directory if it doesn't exist.
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

fwrite(points, file.path(outdir, "HordeumSC_PCoA_points.csv"))
fwrite(segments, file.path(outdir, "HordeumSC_PCoA_segments.csv"))
fwrite(text_dt, file.path(outdir, "HordeumSC_PCoA_text.csv"))

# p <- ggplot(points, aes(x=x, y=y, color=Genotype)) +
#   geom_point(size=1.5, alpha=0.7) +
#   facet_wrap(~Compartment)+
#   geom_segment(data=segments, aes(x = x, y = y, xend = seg_x,
#                                   yend = seg_y, color = Genotype),
#                alpha=0.5) +
#   geom_label(data = text_dt, aes(x = -0.2, y = -0.3, label = text),
#              colour = "black", fill = "grey", alpha = 0.2, size = 8/.pt)+
#   scale_color_manual(values=colors, labels=legend_labels) +
#   guides(color = guide_legend(override.aes = list(linetype = 0))) +
#   labs(
#     x = "PCo 1",
#     y = "PCo 2"
#   ) +
#   main_theme +
#   theme(
#     plot.title = element_text(face="bold", size=8, hjust=0),
#     legend.text = element_markdown(size=8, color="black"),
#     strip.text = element_text(size = 8, colour = "black", face = "bold"),
#     legend.key.size = unit(0.25, "cm")
#   )+
#   NULL
# p
