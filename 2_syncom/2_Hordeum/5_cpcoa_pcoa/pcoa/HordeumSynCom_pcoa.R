# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data.
design <- read.table("../BarleyCSSP_SConly_metadata_NEW.txt", header=T, sep="\t")
asv_table <- read.table("../feature-table_BarleyCSSP_CerealSConly.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, comment.char = "", skip = 1)

# Load required packages.
library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggtext)
library(cowplot)

# Filter design file to only keep genotypes and compartments of interest, and match asv table file accordingly.
design_filtered <- design %>%
  filter(Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2"),
         Compartment %in% c("rhizo","endo")) %>%
  mutate(Compartment = recode(Compartment, "rhizo"="Rhizosphere", "endo"="Root"))

samples_keep <- design_filtered$SampleID
asv_table_filt <- asv_table[, colnames(asv_table) %in% samples_keep]

## Define colours for later graph as well as main theme. Also, indicate that mutant names shall be italic in legend.
colors <- c("WT"="#A9C289","symrk"="#FEDA8B","ccamk"="#FDB366","nsp1"="#C0E4EF","nsp2"="#6EA6CD")

legend_labels <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

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
                    legend.position="right",
                    # legend.background=element_rect(colour="black", fill=NA),
                    legend.background=element_blank(),
                    plot.title = element_text(size=20, hjust=0.9))

# Create function for pcoa analysis and plotting (separately for each compartment).
plot_pcoa_by_compartment <- function(asv_data, design_df){
  compartments <- unique(design_df$Compartment)
  plot_list <- list()
  
  for(comp in compartments){
    # Subset samples
    design_sub <- design_df %>% filter(Compartment == comp)
    asv_sub <- asv_data[, colnames(asv_data) %in% design_sub$SampleID]
    
    # Normalize ASVs
    asv_norm <- apply(asv_sub, 2, function(x) x / sum(x))
    
    # Bray-Curtis & PCoA
    bray <- vegdist(t(asv_norm), method="bray")
    pcoa <- cmdscale(bray, k=2, eig=TRUE)
    points <- as.data.frame(pcoa$points)
    colnames(points) <- c("x","y")
    points <- cbind(points, design_sub[match(rownames(points), design_sub$SampleID), ])
    points$Genotype <- factor(points$Genotype, levels=names(colors))
    
    # Centroids & segments
    centroids <- aggregate(cbind(x,y) ~ Genotype, data=points, FUN=mean)
    segments <- merge(points, setNames(centroids, c("Genotype","seg_x","seg_y")), by="Genotype")
    
    # Plot per compartment
    p <- ggplot(points, aes(x=x, y=y, color=Genotype)) +
      geom_point(size=6, alpha=0.7) +
      geom_segment(data=segments, aes(x=x, y=y, xend=seg_x, yend=seg_y, color=Genotype), alpha=0.5) +
      scale_color_manual(values=colors, labels=legend_labels) +
      guides(color = guide_legend(override.aes = list(size=5))) +
      labs(
        x=paste0("PCoA 1 (", format(100*pcoa$eig[1]/sum(pcoa$eig), digits=4), "%)"),
        y=paste0("PCoA 2 (", format(100*pcoa$eig[2]/sum(pcoa$eig), digits=4), "%)"),
        title=comp
      ) +
      main_theme +
      theme(
        plot.title = element_text(face="bold", size=20, hjust=0),
        legend.text = element_markdown(size=20, color="black")
      )
    
    plot_list[[comp]] <- p
  }
  
  # Combine side-by-side with single legend
  combined <- plot_grid(plot_list[[1]] + theme(legend.position="none"),
                        plot_list[[2]] + theme(legend.position="none"),
                        ncol=2, align="hv")
  legend <- get_legend(plot_list[[1]])
  final <- plot_grid(combined, legend, ncol=2, rel_widths=c(2,0.3))
  
  return(final)
}

# For the pcoa analysis, we want to use two approaches:
# 1. Take all ASVs into account (matched and contaminants).
# 2. Take only matched ASVs into account.

# Run function and make plot for all ASVs.
p_allASVs <- plot_pcoa_by_compartment(asv_table_filt, design_filtered)
ggsave("HordeumSynCom_pcoa_allASVs_byCompartment.pdf", p_allASVs, width=12, height=5)
saveRDS(p_allASVs, "HordeumSynCom_pcoa_allASVs_byCompartment.rds")

# Plot matched ASVs only.
asv_table_matched <- asv_table_filt[grepl("_", rownames(asv_table_filt)), , drop = FALSE]
p_matchedASVs <- plot_pcoa_by_compartment(asv_table_matched, design_filtered)
ggsave("HordeumSynCom_pcoa_matchedASVs_byCompartment.pdf", p_matchedASVs, width=12, height=5)
saveRDS(p_matchedASVs, "HordeumSynCom_pcoa_matchedASVs_byCompartment.rds")


