# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("data.table", "ggplot2", "vegan", "ggtext", "ggpubr","patchwork")

for(pk in pkg){
  library(pk, character.only = T)
}

# Loading data
asv_table <- fread(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv"
)
colnames(asv_table)[1] <- "ASVid"

design <- fread("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt")

# Setting genotype colours
color <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

# Vector for mutants in italic
genotype_labels_legend <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

# Main theme
main_theme <- theme(
  axis.title = element_text(size = 6, colour = "black"),
  axis.text = element_text(size = 6, colour = "black"),
  legend.text = element_markdown(size = 6, colour = "black"),
  legend.title = element_text(size = 6, colour = "black", face = "bold"),
  plot.title = element_text(
    size = 6,
    colour = "black",
    face = "bold",
    hjust = 0.5
  )
)

# CPCoA plots -----------------------------------------------------------------
# Constrained principal coordinate analysis (CPCoA) with plot for each compartment
cond <- c("Root", "Rhizosphere")
for(i in 1:length(cond)){
  # Subsetting design and ASV table
  design_sub <- design[Compartment == cond[i]]
  sub_samples <- design_sub$SampleID
  asv_table_sub <- asv_table[,..sub_samples]
  
  # Relative abundances (RA)
  asv_table_RA <- apply(asv_table_sub, 2, function(x) x/sum(x))
  rownames(asv_table_RA) <- asv_table$ASVid
  
  # Bray-Curtis distances
  bray_curtis <- vegdist(t(asv_table_RA), method = "bray")
  
  # CPCoA
  CPCoA <- capscale(
    t(asv_table_RA) ~ Genotype * Soil,
    data = design_sub,
    add = F,
    sqrt.dist = T,
    distance = "bray"
  )
  
  # Permanova
  set.seed(1762263595)
  nperm <- 999
  perm_test <- anova.cca(CPCoA, permutations = nperm)
  var_expl_tot <- CPCoA$CCA$tot.chi/CPCoA$tot.chi
  p_val <- perm_test["Model", "Pr(>F)"]
  p_val_print <- ifelse(
    p_val == 1 / (nperm + 1),
    paste0("p < ", p_val),
    paste0("p = ", p_val)
  )
  lower_print <- paste0(
    round(var_expl_tot*100, 2), "% of variance; ", p_val_print
  )
  
  CPCo_points <- data.table(
    SampleID = rownames(CPCoA$CCA$wa),
    CPCoA$CCA$wa[, 1:2]
  )
  CPCo_points <- merge(CPCo_points, design_sub, "SampleID")
  var_expl <- CPCoA$CCA$eig/sum(CPCoA$CCA$eig[CPCoA$CCA$eig>0])
  
  # Segments
  CPCo_points[,":="(
    Soil = factor(Soil, levels = c("NPK", "PK", "UF")),
    Genotype = factor(Genotype, levels = names(color))
  )]
  
  # Identifying centroids by soil-compartment combination
  centroids <- CPCo_points[
    ,.(seg_x = mean(CAP1), seg_y = mean(CAP2)), list(Soil, Genotype)
  ]
  segments <- merge(CPCo_points, centroids, by = c("Soil","Genotype"))
  
  # Plots
  CPCoA_plot <- ggplot(
    data = CPCo_points,
    aes(x = CAP1, y = CAP2, colour = Genotype, shape = Soil)
  ) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_segment(
      data = segments,
      alpha = 0.5,
      show.legend = FALSE,
      aes(x = CAP1, y = CAP2, xend = seg_x, yend = seg_y, color = Genotype)
    ) +
    scale_color_manual(values = color, labels = genotype_labels_legend) +
    scale_shape_manual(values = c(15, 17, 16), breaks = c("NPK", "PK", "UF")) +
    theme_bw() +
    main_theme +
    theme(
      legend.key.size = unit(0.25, 'cm'),
      legend.key.spacing.y = unit(0, 'cm')
    ) +
    ggtitle(paste(cond[i], "-", "All soils")) +
    labs(
      x = paste0(
        "CPCo 1 (",
        round(var_expl[1] * 100, 2),
        "%)",
        "\n",
        lower_print
      ),
      y = paste0("CPCo 2 (", round(var_expl[2] * 100, 2), "%)")
    ) +
    NULL
  
  assign( paste("CPCoA_plot", cond[i], sep = "_"), CPCoA_plot )
}

# PCoA plots ------------------------------------------------------------------
# Principal coordinate analysis (PCoA) with plot for each compartment-soil combination
opt <- expand.grid(
  Compartment = c("Root", "Rhizosphere"),
  Soil = c("UF", "PK", "NPK")
)

for(i in 1:nrow(opt)){
  # Subsetting design and ASV table
  current_compartment <- opt$Compartment[i]
  current_soil <- opt$Soil[i]
  design_sub <- design[
    Compartment == current_compartment & Soil == current_soil
  ]
  sub_samples <- design_sub$SampleID
  asv_table_sub <- asv_table[,..sub_samples]
  
  # Relative abundances (RA)
  asv_table_RA <- apply(asv_table_sub, 2, function(x) x/sum(x))
  rownames(asv_table_RA) <- asv_table$ASVid

  # Bray-Curtis distances
  bray_curtis <- vegdist(t(asv_table_RA), method = "bray")
  
  # PCoA
  PCoA <- cmdscale(bray_curtis, k = 2, eig = T)
  PCoA_points <- data.table(SampleID = rownames(PCoA$points), PCoA$points)
  setnames(PCoA_points, c("V1", "V2"), c("PCo1", "PCo2"))
  PCoA_points <- merge(PCoA_points, design_sub, "SampleID")
  var_expl <- PCoA$eig/sum(PCoA$eig[PCoA$eig>0])
  
  PCoA_points[,":="(
    Genotype = factor(Genotype, levels = names(color))
  )]
  
  # Identifying centroids by soil-compartment combination
  centroids <- PCoA_points[,.(seg_x = mean(PCo1), seg_y = mean(PCo2)), Genotype]
  PCoA_points <- merge(PCoA_points, centroids, by = "Genotype")
  
  PCoA_plot <- ggplot(PCoA_points, aes(x = PCo1, y = PCo2, colour = Genotype)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_colour_manual(values = color, labels = genotype_labels_legend) +
    geom_segment(
      aes(xend = seg_x, yend = seg_y),
      alpha = 0.5,
      show.legend = FALSE
    ) +
    ggtitle(paste(current_compartment, "-", current_soil, "soil")) +
    labs(
      x = paste0("PCo 1 (", round(100 * var_expl[1], 2), "%)", "\n", " "),
      y = paste0("PCo 2 (", round(100 * var_expl[2], 2), "%)")
    ) +
    theme_bw() +
    main_theme +
    NULL
  
  assign(
    paste("PCoA", "plot", current_compartment, current_soil, sep = "_"),
    PCoA_plot
  )
}

# Combine plots for the main figure and save.
text <- ggplot() +
  annotate(
    "text",
    x = 1,
    y = 1,
    label = "Hordeum",
    size = 6 / .pt,
    fontface = "bold",
    colour = "black"
  ) +
  theme(
    panel.background = element_rect(fill = "lightgrey", color = NA),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) +
  NULL

CPCoA_PCoA_plot <- text /
  ((CPCoA_plot_Root + theme(legend.position = "none")) +
    (PCoA_plot_Root_UF + theme(legend.position = "none"))) +
  plot_layout(heights = c(0.125, 0.875))

# Saving plot
ggsave(
  filename = "2_figures/HordeumCSSP_AskovSoils_cpcoaRootAll_pcoaRootUF.pdf",
  plot = CPCoA_PCoA_plot,
  width = 9,
  height = 5.5,
  units = "cm"
)

saveRDS(
  object = CPCoA_PCoA_plot,
  file = "1_rds_files/HordeumCSSP_AskovSoils_cpcoaRootAll_pcoaRootUF.rds"
)

# Saving individual plot to use legend later
saveRDS(
  object = CPCoA_plot_Root,
  file = "1_rds_files/HordeumCSSP_AskovSoils_cpcoa_with_legend.rds"
)

# Combine plots for supplementary figures
PCoA_plots <- ggarrange(
  PCoA_plot_Rhizosphere_NPK + theme(legend.position = "none"),
  PCoA_plot_Root_NPK + theme(legend.position = "none"),
  PCoA_plot_Rhizosphere_PK + theme(legend.position = "none"),
  PCoA_plot_Root_PK + theme(legend.position = "none"),
  PCoA_plot_Rhizosphere_UF + theme(legend.position = "none"),
  PCoA_plot_Root_UF + theme(legend.position = "none"),
  ncol = 2, nrow = 3
)

PCoA_plots <- annotate_figure(
  PCoA_plots, top = text_grob("Hordeum", size = 6, face = "bold")
)

# Saving plot
ggsave(
  filename = "2_figures/HordeumCSSP_AskovSoils_pcoa_all.pdf",
  plot = PCoA_plots,
  width = 9,
  height = 15,
  units = "cm"
)

saveRDS(object = PCoA_plots, "1_rds_files/HordeumCSSP_AskovSoils_pcoa_all.rds")

# Combine CPCoA plots for supplementary figures
CPCoA_plots <- ggarrange(
  CPCoA_plot_Rhizosphere + theme(legend.position = "none"),
  CPCoA_plot_Root + theme(legend.position = "none"),
  ncol = 2, nrow = 1
)

CPCoA_plots <- annotate_figure(
  CPCoA_plots, top = text_grob("Hordeum", size = 6, face = "bold")
)

# Saving plots
ggsave(
  filename = "2_figures/HordeumCSSP_AskovSoils_cpcoa_all.pdf",
  plot = CPCoA_plots,
  width = 12,
  height = 5.5,
  units = "cm"
)

saveRDS(
  object = CPCoA_plots,
  file = "1_rds_files/HordeumCSSP_AskovSoils_cpcoa_all.rds"
)
