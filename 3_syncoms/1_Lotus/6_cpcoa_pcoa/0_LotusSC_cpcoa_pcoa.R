# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load plotting functions.
pkg <- c("data.table", "magrittr", "ggplot2", "vegan", "ggtext", "ggpubr", 
         "grid","tidyverse","patchwork")
for(pk in pkg) library(pk, character.only=TRUE)

# Load input files.
asv_table_all <- fread("../1_data/LotusSC_ASVtable_nospike.tsv")
setnames(asv_table_all, "V1", "ASVid")

# design <- fread("../../1_data/LotusSC_metadata.txt")
meta_data <- fread("../1_data/LotusSC_metadata.txt")

# For plotting the beta-diversity based on Bray Curtis dissimilarities, we will use to approaches:
## 1. Take all ASVs into account (matched ASVs and contaminants).
## 2. Take only matched ASVs into account.

# Make filtered asv_table for approach where only matched ASVs are taken into account.
asv_table_matched <- asv_table_all[grepl("Lj", ASVid)]

# Set colours and shapes for plot.
colors <- c("WT"="#A9C289","symrk"="#FEDA8B","ccamk"="#FDB366","nsp1"="#C0E4EF","nsp2"="#6EA6CD")
shapes <- c("Rhizosphere"=15,"Root"=16,"Nodules"=17)

# Make genotype names italic.
genotype_labels_legend <- c("WT"="WT","symrk"="*symrk*","ccamk"="*ccamk*",
                            "nsp1"="*nsp1*","nsp2"="*nsp2*")

# Set plotting parameters.
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 8, color = "black"),
                    # legend.text = element_text(size=8, color = "black"),
                    legend.text = ggtext::element_markdown(size=8, color = "black"),
                    legend.key=element_blank(),
                    axis.title.y = element_text(size = 8),
                    text=element_text(size=8, color="black"),
                    legend.position="right",
                    # legend.background=element_rect(colour="black", fill=NA),
                    legend.background=element_blank(),
                    plot.title = element_text(size=8, hjust=0.9))

# Write the function for making cpcoa and pcoa plots.

run_beta_diversity <- function(asv_table, table_name) {

  ## cpcoa plot
  asv_mat <- as.matrix(asv_table[, -1, with=FALSE])
  rownames(asv_mat) <- asv_table$ASVid
  asv_RA <- apply(asv_mat, 2, function(x) x/sum(x))

  # Bray-Curtis
  bray_curtis <- vegdist(t(asv_RA), method="bray")

  # CPCoA
  CPCoA <- capscale(t(asv_RA) ~ Genotype*Compartment, data=meta_data,
                    add=F, sqrt.dist=T, distance="bray")

  # Permanova
  set.seed(1762263595)
  nperm <- 999
  perm_test <- anova.cca(CPCoA, permutations = nperm)
  var_expl_tot <- CPCoA$CCA$tot.chi/CPCoA$tot.chi
  p_val <- perm_test["Model", "Pr(>F)"]
  p_val_print <- ifelse(p_val == 1/(nperm+1),
                        paste0("p < ", p_val),
                        paste0("p = ", p_val))
  Lower_print <- paste0(round(var_expl_tot*100, 2), "% of variance; ", p_val_print)

  # Sample scores
  CPCo_points <- data.table(SampleID=rownames(CPCoA$CCA$wa), CPCoA$CCA$wa[,1:2])
  CPCo_points <- merge(CPCo_points, meta_data, by="SampleID")
  var_expl <- CPCoA$CCA$eig/sum(CPCoA$CCA$eig[CPCoA$CCA$eig>0])

  CPCo_points[, Genotype := factor(Genotype, levels=names(colors))]
  CPCo_points[, Compartment := factor(Compartment, levels=names(shapes))]

  # Centroids
  centroids <- CPCo_points[, .(seg_x=mean(CAP1), seg_y=mean(CAP2)), by=.(Compartment, Genotype)]
  segments <- merge(CPCo_points, centroids, by=c("Compartment","Genotype"))

  # Plot CPCoA
  CPCoA_plot <- ggplot(CPCo_points, aes(x=CAP1, y=CAP2, colour=Genotype, shape=Compartment)) +
    geom_segment(data=segments, aes(x=CAP1, y=CAP2, xend=seg_x, yend=seg_y, color=Genotype),
                 alpha=0.5, show.legend=FALSE) +
    geom_point(size=5, alpha=0.7) +
    scale_color_manual(values=colors, labels=genotype_labels_legend) +
    scale_shape_manual(values=shapes) +
    labs(x = paste0("CPCo 1 (", round(var_expl[1]*100, 2), "%)", "\n", Lower_print),
         y = paste0("CPCo 2 (", round(var_expl[2]*100, 2), "%)")) +
    main_theme

  # Save CPCoA plot
  ggsave(paste0("LotusSC_cpcoa_", table_name, ".pdf"), CPCoA_plot, width=5, height=5)
  saveRDS(CPCoA_plot, paste0("LotusSC_cpcoa_", table_name, ".rds"))
  saveRDS(CPCoA_plot, paste0("../../3_final_figures/LotusSC_cpcoa_", table_name, ".rds"))

  # pcoa plots (rhizosphere and root).
  for(comp in c("Root","Rhizosphere")) {
    meta_sub <- meta_data[Compartment==comp]
    asv_sub <- asv_table[, c("ASVid", meta_sub$SampleID), with=FALSE]

    rownames(asv_sub) <- asv_sub$ASVid
    asv_sub <- asv_sub[, -1, with=FALSE]
    asv_sub_mat <- as.matrix(asv_sub)
    asv_RA <- apply(asv_sub, 2, function(x) x/sum(x))
    bray_curtis <- vegdist(t(asv_RA), method="bray")

    PCoA <- cmdscale(bray_curtis, k=2, eig=TRUE)
    PCoA_points <- data.table(SampleID = rownames(PCoA$points), PCoA$points)
    setnames(PCoA_points, c("V1", "V2"), c("PCo1", "PCo2"))
    PCoA_points <- merge(PCoA_points, meta_sub, by="SampleID")
    PCoA_points[, Genotype := factor(Genotype, levels=names(colors))]

    centroids <- PCoA_points[, .(seg_x=mean(PCo1), seg_y=mean(PCo2)), by=Genotype]
    PCoA_points <- merge(PCoA_points, centroids, by="Genotype")
    var_expl <- PCoA$eig/sum(PCoA$eig[PCoA$eig>0])

    # Plot PCoA
    p <- ggplot(PCoA_points, aes(x=PCo1, y=PCo2, colour=Genotype)) +
      geom_point(size=3, alpha=0.7) +
      geom_segment(aes(xend=seg_x, yend=seg_y), alpha=0.5, show.legend=FALSE) +
      scale_color_manual(values=colors, labels=genotype_labels_legend) +
      labs(x=paste0("PCo 1 (", round(var_expl[1]*100,2), "%)"),
           y=paste0("PCo 2 (", round(var_expl[2]*100,2), "%)")) +
      ggtitle(comp) +
      main_theme +
      theme(plot.title = element_text(size=8, hjust=0))

    # Save PCoA
    ggsave(paste0("LotusSC_pcoa_", table_name, "_", comp, ".pdf"), p, width=5, height=5)
    saveRDS(p, paste0("LotusSC_pcoa_", table_name, "_", comp, ".rds"))
    saveRDS(p, paste0("../../3_final_figures/LotusSC_pcoa_", table_name, "_", comp, ".rds"))
  }
}

# Run the function for both ASV tables (all ASVs and matched ASVs only).
run_beta_diversity(asv_table_all, "all_ASVs")
run_beta_diversity(asv_table_matched, "matched_ASVs")
