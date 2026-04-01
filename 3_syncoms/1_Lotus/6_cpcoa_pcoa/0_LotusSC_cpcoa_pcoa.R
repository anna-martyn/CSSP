# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("data.table", "ggplot2", "vegan")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Loading data
asv_table_all <- fread("../1_data/LotusSC_ASVtable_nospike.tsv")
setnames(asv_table_all, "V1", "ASVid")

design <- fread("../1_data/LotusSC_metadata.txt")

# Filtered ASV table where only ASVs matched to SynCom sequences are kept
asv_table_matched <- asv_table_all[grepl("Lj", ASVid)]

# Genotype colours
colors <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

# Compartment shapes
shapes <- c("Rhizosphere" = 15, "Root" = 16, "Nodules" = 17)

# Genotype labels to make mutant names italics
genotype_labels_legend <- c(
  "WT" = "WT",
  "symrk" = "*symrk*",
  "ccamk" = "*ccamk*",
  "nsp1" = "*nsp1*",
  "nsp2" = "*nsp2*"
)

# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(size = 6, color = "black"),
  legend.text = ggtext::element_markdown(size = 6, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 6),
  text = element_text(size = 6, color = "black"),
  legend.position = "right",
  legend.background = element_blank(),
  plot.title = element_text(size = 6, hjust = 0.9)
)

# CPCoA and PCoA --------------------------------------------------------------
# Function for CPCoA and PCoA plots
run_beta_diversity <- function(asv_table, table_name) {
  # Relative abundances
  asv_mat <- as.matrix(asv_table[, -1, with = FALSE])
  rownames(asv_mat) <- asv_table$ASVid
  asv_RA <- apply(asv_mat, 2, function(x) x / sum(x))

  # Bray-Curtis
  bray_curtis <- vegdist(t(asv_RA), method = "bray")

  # CPCoA
  cpcoa <- capscale(
    t(asv_RA) ~ Genotype * Compartment,
    data = design,
    add = FALSE,
    sqrt.dist = TRUE,
    distance = "bray"
  )

  # Permanova
  set.seed(1762263595)
  nperm <- 999
  perm_test <- anova.cca(cpcoa, permutations = nperm)
  var_expl_tot <- cpcoa$CCA$tot.chi / cpcoa$tot.chi
  p_val <- perm_test["Model", "Pr(>F)"]
  p_val_print <- ifelse(
    p_val == 1 / (nperm + 1),
    paste0("p < ", p_val),
    paste0("p = ", p_val)
  )
  Lower_print <- paste0(
    round(var_expl_tot * 100, 2),
    "% of variance; ",
    p_val_print
  )

  # Sample scores
  CPCo_points <- data.table(
    SampleID = rownames(cpcoa$CCA$wa),
    cpcoa$CCA$wa[, 1:2]
  )
  CPCo_points <- merge(CPCo_points, design, by = "SampleID")
  var_expl <- cpcoa$CCA$eig / sum(cpcoa$CCA$eig[cpcoa$CCA$eig > 0])

  CPCo_points[, Genotype := factor(Genotype, levels = names(colors))]
  CPCo_points[, Compartment := factor(Compartment, levels = names(shapes))]

  # Centroids
  centroids <- CPCo_points[,
    .(seg_x = mean(CAP1), seg_y = mean(CAP2)),
    by = .(Compartment, Genotype)
  ]
  segments <- merge(CPCo_points, centroids, by = c("Compartment", "Genotype"))

  # Plot CPCoA
  cpcoa_plot <- ggplot(
    CPCo_points,
    aes(x = CAP1, y = CAP2, colour = Genotype, shape = Compartment)
  ) +
    geom_segment(
      data = segments,
      aes(x = CAP1, y = CAP2, xend = seg_x, yend = seg_y, color = Genotype),
      alpha = 0.5,
      show.legend = FALSE
    ) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = colors, labels = genotype_labels_legend) +
    scale_shape_manual(values = shapes) +
    labs(
      x = paste0(
        "CPCo 1 (",
        round(var_expl[1] * 100, 2),
        "%)",
        "\n",
        Lower_print
      ),
      y = paste0("CPCo 2 (", round(var_expl[2] * 100, 2), "%)")
    ) +
    main_theme +
    ggtitle("Lotus") +
    theme(plot.title = element_text(face = "bold", size = 6, hjust = 0))

  # Save CPCoA plot
  ggsave(
    filename = paste0("2_figures/LotusSC_cpcoa_", table_name, ".pdf"),
    plot = cpcoa_plot,
    width = 5,
    height = 5,
    units = "cm"
  )
  saveRDS(
    object = cpcoa_plot,
    file = paste0("1_rds_files/LotusSC_cpcoa_", table_name, ".rds")
  )

  # pcoa plots (rhizosphere and root)
  text_lst <- list(NA, NA)
  pcoa_points_lst <- list(NA, NA)
  names(text_lst) <- c("Root", "Rhizosphere")
  names(pcoa_points_lst) <- c("Root", "Rhizosphere")
  for (comp in c("Root", "Rhizosphere")) {
    # Subsetting compartment
    design_sub <- design[Compartment == comp]
    asv_sub <- asv_table[, c("ASVid", design_sub$SampleID), with = FALSE]

    # Relative abundances
    rownames(asv_sub) <- asv_sub$ASVid
    asv_sub <- asv_sub[, -1, with = FALSE]
    asv_sub_mat <- as.matrix(asv_sub)
    asv_RA <- apply(asv_sub, 2, function(x) x / sum(x))
    
    # Bray-Curtis
    bray_curtis <- vegdist(t(asv_RA), method = "bray")

    # PCoA
    pcoa <- cmdscale(bray_curtis, k = 2, eig = TRUE)
    pcoa_points <- data.table(SampleID = rownames(pcoa$points), pcoa$points)
    setnames(pcoa_points, c("V1", "V2"), c("PCo1", "PCo2"))
    pcoa_points <- merge(pcoa_points, design_sub, by = "SampleID")
    pcoa_points[, Genotype := factor(Genotype, levels = names(colors))]

    # Centroids
    centroids <- pcoa_points[,
      .(seg_x = mean(PCo1), seg_y = mean(PCo2)),
      by = Genotype
    ]
    pcoa_points <- merge(pcoa_points, centroids, by = "Genotype")

    # Variance explained
    var_expl <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])
    text_dt <- data.table(
      Compartment = comp,
      text = paste0(
        round(var_expl[1] * 100, 1),
        "%",
        "-",
        round(var_expl[2] * 100, 1),
        "%"
      )
    )

    # Savings text and points in list
    text_lst[[comp]] <- text_dt
    pcoa_points_lst[[comp]] <- pcoa_points

    # Plot PCoA
    pcoa_plot <- ggplot(
      pcoa_points,
      aes(x = PCo1, y = PCo2, colour = Genotype)
    ) +
      geom_point(size = 3, alpha = 0.7) +
      geom_segment(
        aes(xend = seg_x, yend = seg_y),
        alpha = 0.5,
        show.legend = FALSE
      ) +
      scale_color_manual(values = colors, labels = genotype_labels_legend) +
      labs(
        x = paste0("PCo 1 (", round(var_expl[1] * 100, 2), "%)"),
        y = paste0("PCo 2 (", round(var_expl[2] * 100, 2), "%)")
      ) +
      main_theme +
      ggtitle(comp) +
      theme(plot.title = element_text(size = 6, hjust = 0))

    # Save PCoA
    ggsave(
      filename = paste0(
        "2_figures/LotusSC_pcoa_",
        table_name,
        "_",
        comp,
        ".pdf"
      ),
      plot = pcoa_plot,
      width = 5,
      height = 5,
      units = "cm"
    )
    saveRDS(
      object = pcoa_plot,
      file = paste0("1_rds_files/LotusSC_pcoa_", table_name, "_", comp, ".rds")
    )
  }

  # Combining text and points from both compartments into one data table each
  text_dt <- rbindlist(text_lst)
  pcoa_points <- rbindlist(pcoa_points_lst)
  text_dt[, Host := "Lotus"]
  pcoa_points[, Host := "Lotus"]

  fwrite(
    x = text_dt,
    file = paste0("3_tables/LotusSC_PCoA_text_", table_name, ".csv")
  )
  fwrite(
    x = pcoa_points,
    file = paste0("3_tables/LotusSC_PCoA_points_segments_", table_name, ".csv")
  )
}

# Run function (all ASVs and matched ASVs)
run_beta_diversity(asv_table = asv_table_all, table_name = "all_ASVs")
run_beta_diversity(asv_table = asv_table_matched, table_name = "matched_ASVs")
