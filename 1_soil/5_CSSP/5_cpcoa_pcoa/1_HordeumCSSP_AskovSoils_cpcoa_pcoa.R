# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
pkg <- c("data.table", "magrittr", "ggplot2", "vegan", "ggtext", "ggpubr", 
         "grid","tidyverse","patchwork")

for(pk in pkg){
  library(pk, character.only = T)
}

# Load the input files.
asv_table <- read.table(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  skip = 1,
  comment.char = ""
)
setDT(asv_table) 

meta_data <- fread("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt")

# Set the colours for all genotypes.
color <- c("WT" = "#A9C289", "symrk" = "#FEDA8B", "ccamk" = "#FDB366",
           "nsp1" = "#C0E4EF", "nsp2" = "#6EA6CD")

# Make the mutant labels italic.
genotype_labels_legend <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

# Set the main theme for the plots.
main_theme <- theme(
  axis.title = element_text(size = 6, colour = "black"),
  axis.text = element_text(size = 6, colour = "black"),
  legend.text = element_markdown(size = 6, colour = "black"),
  legend.title = element_text(size = 6, colour = "black", face = "bold"),
  plot.title = element_text(size = 6, colour = "black", face = "bold", 
                            hjust = 0.5)
)

# Perform constrained principal coordinate analysis (CPCoA) for each compartment, and make the corresponding plots.
cond <- c("Root", "Rhizosphere")
for(i in 1:length(cond)){
  meta_data_sub <- meta_data[Compartment == cond[i]]
  sub_samples <- meta_data_sub$SampleID
  asv_table_sub <- asv_table[,..sub_samples]
  
  # Convert reads to relative abundances for each ASV.
  asv_table_RA <- apply(asv_table_sub, 2, function(x) x/sum(x))
  rownames(asv_table_RA) <- asv_table$ASVid
  
  # Calculate Bray-Curtis distances.
  bray_curtis <- vegdist(t(asv_table_RA), method = "bray")
  
  CPCoA <- capscale(t(asv_table_RA) ~ Genotype*Soil, data = meta_data_sub,
                    add = F, sqrt.dist = T, distance = "bray")
  
  # Permanova.
  set.seed(1762263595)
  nperm <- 999
  perm_test <- anova.cca(CPCoA, permutations = nperm)
  var_expl_tot <- CPCoA$CCA$tot.chi/CPCoA$tot.chi
  p_val <- perm_test["Model", "Pr(>F)"]
  p_val_print <- ifelse(p_val == 1/(nperm+1),
                        paste0("p < ", p_val),
                        paste0("p = ", p_val))
  Lower_print <- paste0(
    round(var_expl_tot*100, 2), "% of variance; ", p_val_print
  )
  
  CPCo_points <- data.table(SampleID = rownames(CPCoA$CCA$wa), 
                            CPCoA$CCA$wa[,1:2])
  CPCo_points <- merge(CPCo_points, meta_data_sub, "SampleID")
  var_expl <- CPCoA$CCA$eig/sum(CPCoA$CCA$eig[CPCoA$CCA$eig>0])
  
  # Segments.
  CPCo_points[,":="(
    Soil = factor(Soil, levels = c("NPK", "PK", "UF")),
    Genotype = factor(Genotype, levels = names(color))
  )]
  
  # Calculate centroids per soil-compartment combination.
  centroids <- CPCo_points[
    ,.(seg_x = mean(CAP1), seg_y = mean(CAP2)), list(Soil, Genotype)
  ]
  segments <- merge(CPCo_points, centroids, by = c("Soil","Genotype"))
  
  # Make the plots.
  ggplot(data = CPCo_points, 
         aes(x = CAP1, y = CAP2, colour = Genotype, shape = Soil))+
    geom_point(size = 1.5, alpha = 0.7)+
    geom_segment(
      data = segments, alpha = 0.5, show.legend = FALSE,
      aes(x = CAP1, y = CAP2, xend = seg_x, yend = seg_y, color = Genotype)
    ) +
    scale_color_manual(values = color, labels = genotype_labels_legend)+
    scale_shape_manual(values = c(15, 17, 16), breaks = c("NPK", "PK", "UF"))+
    theme_bw()+
    main_theme+
    theme(legend.key.size = unit(0.25, 'cm'),
          legend.key.spacing.y = unit(0, 'cm'))+
    ggtitle(paste(cond[i], "-", "All soils"))+
    labs(
      x = paste0("CPCo 1 (", round(var_expl[1]*100, 2), "%)", "\n", Lower_print),
      y = paste0("CPCo 2 (", round(var_expl[2]*100, 2), "%)")
    )+
    NULL -> g
  
  assign( paste("CPCoA_plot", cond[i], sep = "_"), g )
}

# Perform principal coordinate analysis (PCoA) for ech compartment-soil combination, and make the corresponding plots.
Opt <- expand.grid(Compartment = c("Root", "Rhizosphere"),
                   Soil = c("UF", "PK", "NPK"))

for(i in 1:nrow(Opt)){
  meta_data_sub <- meta_data[
    Compartment == Opt$Compartment[i] & Soil == Opt$Soil[i]
  ]
  sub_samples <- meta_data_sub$SampleID
  asv_table_sub <- asv_table[,..sub_samples]
  
  asv_table_RA <- apply(asv_table_sub, 2, function(x) x/sum(x))
  rownames(asv_table_RA) <- asv_table$ASVid
  bray_curtis <- vegdist(t(asv_table_RA), method = "bray")
  
  PCoA <- cmdscale(bray_curtis, k = 2, eig = T)
  PCoA_points <- data.table(SampleID = rownames(PCoA$points), PCoA$points)
  setnames(PCoA_points, c("V1", "V2"), c("PCo1", "PCo2"))
  PCoA_points <- merge(PCoA_points, meta_data_sub, "SampleID")
  var_expl <- PCoA$eig/sum(PCoA$eig[PCoA$eig>0])
  
  PCoA_points[,":="(
    Genotype = factor(Genotype, levels = names(color))
  )]
  # Calculate centroids per soil-compartment combination.
  centroids <- PCoA_points[,.(seg_x = mean(PCo1), seg_y = mean(PCo2)), Genotype]
  PCoA_points <- merge(PCoA_points, centroids, by = "Genotype")
  
  ggplot(PCoA_points, aes(x = PCo1, y = PCo2, colour = Genotype))+
    geom_point(alpha = 0.7, size = 1.5) + 
    scale_colour_manual(values = color, labels = genotype_labels_legend) +
    geom_segment(
      aes(xend = seg_x, yend = seg_y), alpha = 0.5, show.legend = FALSE
    ) +
    ggtitle( paste(Opt$Compartment[i], "-", Opt$Soil[i], "soil") ) +
    labs(
      x = paste0("PCo 1 (", round(100*var_expl[1], 2), "%)", "\n", " "),
      y = paste0("PCo 2 (", round(100*var_expl[2], 2), "%)")
    ) +
    theme_bw()+
    main_theme +
    NULL -> g
  
  assign(paste("PCoA", "plot", Opt$Compartment[i], Opt$Soil[i], sep = "_"), g)
}

# Combine plots for the main figure and save.
text <- ggplot() +
  annotate("text", x = 1, y = 1, label = "Hordeum", size = 6/.pt,
           fontface = "bold", colour = "black") +
  theme(
    panel.background = element_rect(fill = "lightgrey", color = NA),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) +
  NULL

text/((CPCoA_plot_Root + theme(legend.position = "none"))+
        (PCoA_plot_Root_UF + theme(legend.position = "none")))+
  plot_layout(heights = c(0.125,0.875)) -> CPCoA_PCoA_plot

ggsave("HordeumCSSP_AskovSoils_cpcoaRootAll_pcoaRootUF.pdf", CPCoA_PCoA_plot, width = 9, height = 5.5, units = "cm")
saveRDS(CPCoA_PCoA_plot, "HordeumCSSP_AskovSoils_cpcoaRootAll_pcoaRootUF.rds")
saveRDS(CPCoA_PCoA_plot, "../8_final_figures/HordeumCSSP_AskovSoils_cpcoaRootAll_pcoaRootUF.rds")

## Also save an inidivudal plot as we'll use the legend of the plot for the final figure later.
saveRDS(CPCoA_plot_Root, "HordeumCSSP_AskovSoils_cpcoa_with_legend.rds")
saveRDS(CPCoA_plot_Root, "../8_final_figures/HordeumCSSP_AskovSoils_cpcoa_with_legend.rds")

# Combine plots for supplementary figures and save.
ggarrange(
  PCoA_plot_Rhizosphere_NPK + theme(legend.position = "none"),
  PCoA_plot_Root_NPK + theme(legend.position = "none"),
  PCoA_plot_Rhizosphere_PK + theme(legend.position = "none"),
  PCoA_plot_Root_PK + theme(legend.position = "none"),
  PCoA_plot_Rhizosphere_UF + theme(legend.position = "none"),
  PCoA_plot_Root_UF + theme(legend.position = "none"),
  ncol = 2, nrow = 3
) -> PCoA_plots

PCoA_plots <- annotate_figure(
  PCoA_plots, top = text_grob("Hordeum", size = 6, face = "bold")
)

ggsave("HordeumCSSP_AskovSoils_pcoa_all.pdf", PCoA_plots, width = 9, height = 15, units = "cm")
saveRDS(PCoA_plots, "HordeumCSSP_AskovSoils_pcoa_all.rds")
saveRDS(PCoA_plots, "../8_final_figures/HordeumCSSP_AskovSoils_pcoa_all.rds")

ggarrange(
  CPCoA_plot_Rhizosphere + theme(legend.position = "none"),
  CPCoA_plot_Root + theme(legend.position = "none"),
  ncol = 2, nrow = 1
) -> CPCoA_plots

CPCoA_plots <- annotate_figure(
  CPCoA_plots, top = text_grob("Hordeum", size = 6, face = "bold")
)

ggsave("HordeumCSSP_AskovSoils_cpcoa_all.pdf", CPCoA_plots, width = 12, height = 5.5, units = "cm")
saveRDS(CPCoA_plots, "HordeumCSSP_AskovSoils_cpcoa_all.rds")
saveRDS(CPCoA_plots, "../8_final_figures/HordeumCSSP_AskovSoils_cpcoa_all.rds")
