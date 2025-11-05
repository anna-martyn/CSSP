# Load packages and export data ------------------------------------------------
pkg <- c("data.table", "magrittr", "ggplot2", "vegan", "ggtext", "ggpubr", 
         "grid")

for(pk in pkg){
  library(pk, character.only = T)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ASV_table <- fread(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv"
)
setnames(ASV_table, "V1", "ASV_ID")
meta_data <- fread(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", drop = c(5,7,8)
)
taxonomy <- fread("../../1_data/1_Lotus/taxonomy_clean.tsv")
setnames(taxonomy, "Feature", "ASV_ID")

color <- c("WT" = "#A9C289", "symrk" = "#FEDA8B", "ccamk" = "#FDB366",
           "nsp1" = "#C0E4EF", "nsp2" = "#6EA6CD")

genotype_labels_legend <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

# Main theme -------------------------------------------------------------------
main_theme <- theme(
  axis.title = element_text(size = 8, colour = "black"),
  axis.text = element_text(size = 8, colour = "black"),
  legend.text = element_markdown(size = 8, colour = "black"),
  legend.title = element_text(size = 8, colour = "black", face = "bold"),
  plot.title = element_text(size = 8, colour = "black", face = "bold", 
                            hjust = 0.5),
)

# CPCoA ------------------------------------------------------------------------
cond <- c("Root", "Rhizosphere")
for(i in 1:length(cond)){
  meta_data_sub <- meta_data[Compartment == cond[i]]
  sub_samples <- meta_data_sub$SampleID
  ASV_table_sub <- ASV_table[,..sub_samples]
  
  ASV_table_RA <- apply(ASV_table_sub, 2, function(x) x/sum(x))
  rownames(ASV_table_RA) <- ASV_table$ASV_ID
  bray_curtis <- vegdist(t(ASV_table_RA), method = "bray")
  
  CPCoA <- capscale(t(ASV_table_RA) ~ Genotype*Soil, data = meta_data_sub,
                    add = F, sqrt.dist = T, distance = "bray")
  
  # Permanova
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
  
  # Segments
  CPCo_points[,":="(
    Soil = factor(Soil, levels = c("UF", "PK", "NPK")),
    Genotype = factor(Genotype, levels = names(color))
  )]
  
  # Calculate centroids per Soil x Compartment.
  centroids <- CPCo_points[
    ,.(seg_x = mean(CAP1), seg_y = mean(CAP2)), list(Soil, Genotype)
  ]
  segments <- merge(CPCo_points, centroids, by = c("Soil","Genotype"))
  
  ggplot(data = CPCo_points, 
         aes(x = CAP1, y = CAP2, colour = Genotype, shape = Soil))+
    geom_point(size = 1.5, alpha = 0.7)+
    geom_segment(
      data = segments, alpha = 0.5, show.legend = FALSE,
      aes(x = CAP1, y = CAP2, xend = seg_x, yend = seg_y, color = Genotype)
    ) +
    scale_color_manual(values = color, labels = genotype_labels_legend)+
    theme_bw()+
    main_theme+
    ggtitle(paste(cond[i], "-", "All soils"))+
    labs(
      x = paste0("CPCo 1 (", round(var_expl[1]*100, 2), "%)", "\n", Lower_print),
      y = paste0("CPCo 2 (", round(var_expl[2]*100, 2), "%)")
    )+
    NULL -> g
  
  assign( paste("CPCoA_plot", cond[i], sep = "_"), g )
}

# PCoA -------------------------------------------------------------------------
Opt <- expand.grid(Compartment = c("Root", "Rhizosphere"),
                   Soil = c("UF", "PK", "NPK"))

for(i in 1:nrow(Opt)){
  meta_data_sub <- meta_data[
    Compartment == Opt$Compartment[i] & Soil == Opt$Soil[i]
  ]
  sub_samples <- meta_data_sub$SampleID
  ASV_table_sub <- ASV_table[,..sub_samples]
  
  ASV_table_RA <- apply(ASV_table_sub, 2, function(x) x/sum(x))
  rownames(ASV_table_RA) <- ASV_table$ASV_ID
  bray_curtis <- vegdist(t(ASV_table_RA), method = "bray")
  
  PCoA <- cmdscale(bray_curtis, k = 2, eig = T)
  PCoA_points <- data.table(SampleID = rownames(PCoA$points), PCoA$points)
  setnames(PCoA_points, c("V1", "V2"), c("PCo1", "PCo2"))
  PCoA_points <- merge(PCoA_points, meta_data_sub, "SampleID")
  var_expl <- PCoA$eig/sum(PCoA$eig[PCoA$eig>0])
  
  PCoA_points[,":="(
    Genotype = factor(Genotype, levels = names(color))
  )]
  # Calculate centroids per Soil x Compartment.
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

# Main figure ------------------------------------------------------------------
text <- ggplot() +
  annotate("text", x = 1, y = 1, label = "Lotus", size = 8/.pt,
           fontface = "bold", colour = "black") +
  theme(
    # plot.background = element_rect(fill = "grey", color = NA),
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
  
saveRDS(CPCoA_PCoA_plot, "CPCoA_PCoA_plot_Lj.rds")

# Supplementary figures --------------------------------------------------------
ggarrange(
  PCoA_plot_Root_UF + theme(legend.position = "none"),
  PCoA_plot_Rhizosphere_UF + theme(legend.position = "none"),
  PCoA_plot_Root_PK + theme(legend.position = "none"),
  PCoA_plot_Rhizosphere_PK + theme(legend.position = "none"),
  PCoA_plot_Root_NPK + theme(legend.position = "none"),
  PCoA_plot_Rhizosphere_NPK + theme(legend.position = "none"),
  ncol = 2, nrow = 3
) -> PCoA_plots

PCoA_plots <- annotate_figure(
  PCoA_plots, top = text_grob("Lotus", size = 8, face = "bold")
)

saveRDS(PCoA_plots, "PCoA_plots_Lj.rds")
