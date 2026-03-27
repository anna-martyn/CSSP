## Original script by Ib Thorsgaard Jensen (Aalborg University),
# modified and extended by Anna Martyn (amartyn@mpipz.mpg.de)

# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Loading packages
pkg <- c(
  "ggplot2",
  "tibble",
  "Maaslin2",
  "ggh4x",
  "dplyr",
  "tidyr",
  "patchwork"
)
for(pk in pkg) library(pk, character.only = TRUE)

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading file with structural zero function
source("Structural_zeros.R")

# Loading data
design <- read.table(
  "../1_data/LotusSC_metadata.txt",
  header = TRUE,
  sep = "\t"
)
asv_table <- read.table(
  "../1_data/LotusSC_ASVtable.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1
)
taxonomy <- read.table(
  "../1_data/LjSC_taxonomy.txt",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)

# Keep only ASVs matched to SynCom
asv_table_matched <- asv_table[grepl("Lj", rownames(asv_table)), ]

# Library sizes (needed for Structural_zeros3)
design$library_size <- colSums(asv_table_matched[, design$SampleID])

# Splitting data by compartment
samples_rhizo <- design$SampleID[design$Compartment == "Rhizosphere"]
samples_root  <- design$SampleID[design$Compartment == "Root"]

asv_table_rhizo <- asv_table_matched[, samples_rhizo, drop = FALSE]
asv_table_root <- asv_table_matched[, samples_root, drop = FALSE]

design_rhizo <- design %>%
  filter(SampleID %in% samples_rhizo) %>%
  column_to_rownames("SampleID")

design_root <- design %>%
  filter(SampleID %in% samples_root) %>%
  column_to_rownames("SampleID")

# Genotype factor levels
design_rhizo$Genotype <- factor(
  design_rhizo$Genotype, levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)
design_root$Genotype  <- factor(
  design_root$Genotype, levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)

# Initialising results table
results_rhizo <- data.frame(matrix(
  NA,
  nrow = nrow(asv_table_matched),
  ncol = 8
))
rownames(results_rhizo) <- rownames(asv_table_matched)
colnames(results_rhizo) <- c(
  paste("Lfc", c("symrk", "ccamk", "nsp1", "nsp2"), sep = "_"),
  paste("DA", c("symrk", "ccamk", "nsp1", "nsp2"), sep = "_")
)
results_root <- results_rhizo

# Structural zero analysis -----------------------------------------------------
da_sz_rhizo <- Structural_zeros3(
  asv_table_rhizo,
  design_rhizo,
  group = "Genotype",
  ref = "WT",
  min_reads = 20,
  min_present_reps = 2
)

da_sz_root <- Structural_zeros3(
  asv_table_root,
  design_root,
  group = "Genotype",
  ref = "WT",
  min_reads = 20,
  min_present_reps = 2
)

# Differential abundance analysis using Maaslin2 -------------------------------
# Rhizosphere
da_maa_rhizo <- Maaslin2(
  input_data = asv_table_rhizo,
  input_metadata = design_rhizo,
  output = "Maaslin2_rhizo",
  fixed_effects = "Genotype",
  min_prevalence = 0.1,
  plot_heatmap = F,
  plot_scatter = F
)
res <- da_maa_rhizo$results

for(g in c("symrk", "ccamk", "nsp1", "nsp2")){
  res_g <- res[res$value == g, ]
  cols <- paste(c("Lfc", "DA"), g, sep = "_")
  rows <- res_g$feature
  results_rhizo[rows, cols[1]] <- res_g$coef
  results_rhizo[rows, cols[2]] <- (res_g$qval < 0.05) * sign(res_g$coef)
}

# Updating results with structural zeros
all_da_rhizo <- Reduce("union", da_sz_rhizo$struc_zero_DA)
struc_zero_wt <- da_sz_rhizo$struc_zero_table[all_da_rhizo, 1]
struc_zero_mutants <- da_sz_rhizo$struc_zero_table[all_da_rhizo, -1]
da_sz_sign_rhizo <- struc_zero_wt - struc_zero_mutants

for(g in c("symrk", "ccamk", "nsp1", "nsp2")){
  idx <- da_sz_rhizo$struc_zero_DA[[g]]
  results_rhizo[idx, paste0("DA_", g)] <- da_sz_sign_rhizo[idx, g]
  results_rhizo[idx, paste0("Lfc_", g)] <- da_sz_sign_rhizo[idx, g] * Inf
}

# Root
da_maa_root <- Maaslin2(
  input_data = asv_table_root,
  input_metadata = design_root,
  output = "Maaslin2_root",
  fixed_effects = "Genotype",
  min_prevalence = 0.1,
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

res <- da_maa_root$results
res %>%
  group_by(feature) %>%
  summarise(DA = any(qval < 0.05)) %>%
  print(n = Inf)

for(g in c("symrk", "ccamk", "nsp1", "nsp2")){
  res_g <- res[res$value == g,]
  cols <- paste(c("Lfc", "DA"), g, sep = "_")
  rows <- res_g$feature
  results_root[rows, cols[1]] <- res_g$coef
  results_root[rows, cols[2]] <- (res_g$qval < 0.05) * sign(res_g$coef)
}

# Updating results with structural zeros
all_da_root <- Reduce("union", da_sz_root$struc_zero_DA)
struc_zero_wt <- da_sz_root$struc_zero_table[all_da_root, 1]
struc_zero_mutants <- da_sz_root$struc_zero_table[all_da_root, -1]
da_sz_sign_root <- struc_zero_wt - struc_zero_mutants

for(g in c("symrk", "ccamk", "nsp1", "nsp2")){
  idx <- da_sz_root$struc_zero_DA[[g]]
  results_root[idx, paste0("DA_", g)] <- da_sz_sign_root[idx, g]
  results_root[idx, paste0("Lfc_", g)] <- da_sz_sign_root[idx, g] * Inf
}

# Merging results with taxonomy and relative abundance info
RA_rhizo <- t(t(asv_table_rhizo)/colSums(asv_table_rhizo))
RA_root  <- t(t(asv_table_root)/colSums(asv_table_root))

results_rhizo <- merge(
  results_rhizo,
  taxonomy,
  by.x = 0,
  by.y = "ASVid"
)
colnames(results_rhizo)[1] <- "ASVid"

results_rhizo <- merge(
  x = results_rhizo,
  y = RA_rhizo,
  by.x = "ASVid",
  by.y = 0
)

results_root <- merge(
  x = results_root,
  y = taxonomy,
  by.x = 0,
  by.y = "ASVid"
)
colnames(results_root)[1] <- "ASVid"

results_root <- merge(
  x = results_root,
  y = RA_root,
  by.x = "ASVid",
  by.y = 0
)

# ASVs to keep in visualisation
rhizo_any <- apply(results_rhizo[,6:9] != 0, 1, any)
root_any <- apply(results_root[,6:9] != 0, 1, any)
rhizo_any[is.na(rhizo_any)] <- FALSE
root_any[is.na(root_any)] <- FALSE
isolate_keep <- results_rhizo$ASVid[rhizo_any|root_any]

# Saving results
write.csv(results_rhizo, file = "DA_SynCom_Lotus_rhizo.csv")
write.csv(results_root, file = "DA_SynCom_Lotus_rhizo.csv")

# Plot highlighting DA ASVs ----------------------------------------------------

# Mean RA for WT in Rhizosphere and Root
asv_table_RA <- sweep(asv_table_matched, 2, colSums(asv_table_matched), "/")
asv_RA_long <- as.data.frame(asv_table_RA) %>%
  rownames_to_column("ASVid") %>%
  pivot_longer(cols=-ASVid, names_to="SampleID", values_to="RA") %>%
  left_join(design %>% select(SampleID, Compartment, Genotype), by="SampleID")

# Filtering for WT
asv_RA_wt <- asv_RA_long %>%
  filter(Genotype=="WT") %>%
  group_by(ASVid, Compartment) %>%
  summarise(mean_RA = mean(RA, na.rm=TRUE), .groups="drop")

# Adding taxonomy
asv_RA_wt <- asv_RA_wt %>%
  left_join(taxonomy %>% select(ASVid, order), by="ASVid") %>%
  mutate(order = ifelse(is.na(order), "Unknown", order))

# Ordering ASVs by taxonomic order
asv_order_levels <- asv_RA_wt %>%
  distinct(ASVid, order) %>%
  arrange(order) %>%
  pull(ASVid)

asv_RA_wt$ASVid <- factor(asv_RA_wt$ASVid, levels = asv_order_levels)

# Compartment factor levels
asv_RA_wt$Compartment <- factor(
  asv_RA_wt$Compartment,
  levels = c("Rhizosphere", "Root", "Nodules")
)

## Taxonomic color bar ---------------------------------------------------------
colors <- read.table(
  "../../../0_files/Bacterial_order_colors.csv",
  header = TRUE,
  sep = ",",
  comment.char = ""
)

tax_bar <- asv_RA_wt %>%
  distinct(ASVid, order) %>% 
  filter(ASVid %in% isolate_keep)

p_tax <- ggplot(tax_bar, aes(x = ASVid, y = 1, fill = order)) +
  geom_tile() +
  scale_fill_manual(values = colors$Color, breaks = colors$Order) +
  theme_void() +
  labs(fill = "Bacterial order") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(color = "black", size = 6),
    legend.title = element_text(color = "black", size = 6),
    legend.key.size = unit(0.25, 'cm'),
    legend.key.spacing.y = unit(0, 'cm'),
    plot.margin = margin(t = 0.25, r = 0, b = 0.25, l = 0, unit = "lines")
    # plot.margin = margin(c(0.25, 0, 0.25, 0), unit = "lines")
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  NULL

## Bubble plot preparation -----------------------------------------------------
# Long form with selected columns
da_rhizo <- results_rhizo %>%
  select(ASVid, DA_symrk, DA_ccamk, DA_nsp1, DA_nsp2) %>%
  pivot_longer(cols = -ASVid, names_to = "Genotype", values_to = "DA") %>%
  mutate(
    Compartment = "Rhizosphere",
    Genotype = recode(
      Genotype,
      DA_symrk = "symrk",
      DA_ccamk = "ccamk",
      DA_nsp1 = "nsp1",
      DA_nsp2 = "nsp2"
    )
  )

da_root <- results_root %>%
  select(ASVid, DA_symrk, DA_ccamk, DA_nsp1, DA_nsp2) %>%
  pivot_longer(cols = -ASVid, names_to = "Genotype", values_to = "DA") %>%
  mutate(
    Compartment = "Root",
    Genotype = recode(
      Genotype,
      DA_symrk = "symrk",
      DA_ccamk = "ccamk",
      DA_nsp1 = "nsp1",
      DA_nsp2 = "nsp2"
    )
  )

# Combining rhizosphere and root tables
da_all <- bind_rows(da_rhizo, da_root)

# Setting factor levels
da_all$ASVid <- factor(da_all$ASVid, levels = asv_order_levels)
da_all$Genotype <- factor(
  da_all$Genotype,
  levels = c("symrk", "ccamk", "nsp1", "nsp2")
)
da_all$Genotype <- factor(
  da_all$Genotype,
  levels = rev(c("symrk", "ccamk", "nsp1", "nsp2"))
)

# Replacing NA with 0 to get blank entries in bubble plot
da_all <- da_all %>% 
  replace_na(list(DA = 0))

da_all$DA <- factor(da_all$DA, levels = c(-1, 0, 1))
da_colors <- c("-1" = "darkblue", "0" = "white", "1" = "red")

# Changing names for appropriate line breaks in plot
da_all$Compartment[da_all$Compartment == "Rhizosphere"] <- "Rhizo-\nsphere"
da_all$Compartment <- factor(
  da_all$Compartment,
  levels = c("Rhizo-\nsphere", "Root")
)

# Keep only ASVs that are DA in at least on mutant for display
da_all <- da_all %>% filter(ASVid %in% isolate_keep)

p_bubble <- ggplot(da_all, aes(x = ASVid, y = Genotype, fill = DA)) +
  geom_point(shape = 21, size = 2, color = "black") +
  scale_fill_manual(
    values = da_colors,
    labels = c("-1" = "Depleted", "0" = "Non-significant", "1" = "Enriched")
  ) +
  facet_grid(Compartment ~ ., scales = "free_x", switch = "y") +
  labs(fill = "Relative abundance mutant vs. WT") +
  labs(y = "Differencial abundance\nin mutants") +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6,
      color = "black"
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 6),
    strip.placement = "outside",
    # strip.background = element_rect(fill = "grey80", color = "grey50"),
    strip.text.y.left = element_text(
      color = "black",
      size = 6,
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.y = element_text(color = "black", size = 6),
    legend.text = element_text(color = "black", size = 6),
    legend.title = element_text(color = "black", size = 6),
    legend.position = "bottom",
    strip.background = element_rect(colour = NA),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "lines"),
    panel.spacing = unit(0.4, "lines")
  ) +
  scale_y_discrete(
    labels = c(
      "symrk" = expression(italic("symrk")),
      "ccamk" = expression(italic("ccamk")),
      "nsp1" = expression(italic("nsp1")),
      "nsp2" = expression(italic("nsp2"))
    )
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  NULL

## Mean RA bar plots -----------------------------------------------------------
### All compartments -----------------------------------------------------------

# Changing names for appropriate line breaks in plot
asv_RA_wt$Compartment <- as.character(asv_RA_wt$Compartment)
asv_RA_wt$Compartment[
  asv_RA_wt$Compartment == "Rhizosphere"
] <- "Rhizo-\nsphere"
asv_RA_wt$Compartment[asv_RA_wt$Compartment == "Nodules"] <- "Nod-\nules"
asv_RA_wt$Compartment <- factor(
  asv_RA_wt$Compartment,
  levels = c("Rhizo-\nsphere", "Root", "Nod-\nules")
)

# Keep only ASVs that are DA in at least on mutant for display
asv_RA_wt <- asv_RA_wt %>% filter(ASVid %in% isolate_keep)

p_RA <- ggplot(asv_RA_wt, aes(x = ASVid, y = mean_RA)) +
  geom_bar(stat = "identity", fill = "grey50") +
  facet_wrap(
    ~Compartment,
    ncol = 1,
    scales = "free_y",
    strip.position = "left",
    space = "free_y"
  ) +
  labs(y = "Mean relative\nabundance in WT") +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ggtitle("Lotus") +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(color = "black", size = 6, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", size = 6),
    strip.text = element_text(color = "black", size = 6, face = "bold"),
    strip.placement = "outside",
    strip.background = element_rect(colour = NA),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "lines"),
    panel.spacing = unit(0.4, "lines")
  ) +
  force_panelsizes(cols = c(1, 1, 1), rows = c(1, 1, 0.7)) +
  # facetted_pos_scales(
  #   y = list(
  #     Compartment == "Nod-\nules" ~ scale_y_continuous(breaks = c(0, 0.25, 0.5))
  #   )
  # )+
  facetted_pos_scales(
    y = list(
      Compartment == "Nod-\nules" ~ scale_y_continuous(
        breaks = c(0, 0.25, 0.5),
        limits = c(0, 0.6),
        expand = c(0, 0)
      ),
      Compartment == "Rhizo-\nsphere" ~ scale_y_continuous(
        limits = c(0, 0.33),
        expand = c(0, 0)
      ),
      Compartment == "Root" ~ scale_y_continuous(
        limits = c(0, 0.33),
        expand = c(0, 0)
      )
    )
  ) +
  NULL

### No nodules -----------------------------------------------------------------

# Filter out Nodules and rename Rhizo-\nsphere back to Rhizosphere again.
# asv_RA_WT_no_nod <- asv_RA_WT %>%
#   filter(Compartment != "Nod-\nules") %>%
#   mutate(Compartment = recode(Compartment, "Rhizo-\nsphere" = "Rhizosphere"),
#          Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")))

# Removing nodule samples and setting factor levels
asv_RA_wt_no_nod <- asv_RA_wt %>%
  filter(Compartment != "Nod-\nules") %>%
  mutate(
    Compartment = factor(Compartment, levels = c("Rhizo-\nsphere", "Root"))
  )

p_RA_no_nod <- ggplot(asv_RA_wt_no_nod, aes(x = ASVid, y = mean_RA)) +
  geom_bar(stat = "identity", fill = "grey50") +
  facet_wrap(
    ~Compartment,
    ncol = 1,
    scales = "free_y",
    strip.position = "left",
    space = "free_y"
  ) +
  labs(y = "Mean relative\nabundance in WT") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ggtitle("Lotus") +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(color = "black", size = 6, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(color = "black", size = 6),
    strip.text = element_text(color = "black", size = 6, face = "bold"),
    strip.placement = "outside",
    strip.background = element_rect(colour = NA),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "lines"),
    panel.spacing = unit(0.4, "lines")
  ) +
  force_panelsizes(cols = c(1, 1), rows = c(1, 1)) +
  facetted_pos_scales(
    y = list(
      Compartment == "Rhizo-\nsphere" ~ scale_y_continuous(
        limits = c(0, 0.33),
        expand = c(0, 0)
      ),
      Compartment == "Root" ~ scale_y_continuous(
        limits = c(0, 0.33),
        expand = c(0, 0)
      )
    )
  )

## Combining plots -------------------------------------------------------------
# once with all compartments, once with rhizo + root only


# Removing individual legends from plots
p_tax_clean <- p_tax + theme(legend.position = "none")
# p_bubble_clean <- p_bubble + theme(legend.position = "none")

# Combining plots vertically
main_plot_all <- p_RA /
  p_tax_clean /
  (p_bubble + theme(legend.position = "none")) +
  plot_layout(heights = c(0.3, 0.05, 0.65))

main_plot_no_nod <- p_RA_no_nod /
  p_tax_clean /
  (p_bubble + theme(legend.position = "none")) +
  plot_layout(heights = c(0.35, 0.05, 0.6))

# lgd_bubble <- ggpubr::get_legend(p_bubble, position = "bottom")
# lgd_tax <- ggpubr::get_legend(p_tax, position = "bottom")
# lgd <- plot_grid(lgd_bubble, lgd_tax, ncol = 2)
#
# plot_grid(final_plot, lgd, rel_heights = c(0.9, 0.1), ncol = 1)

# Saving plots
ggsave(
  "LotusSynCom_DA_withNodule.pdf",
  plot = main_plot_all,
  width = 21,
  height = 20,
  units = "cm"
)
saveRDS(main_plot_all, file = "LotusSynCom_DA_withNodule.rds")
saveRDS(
  main_plot_all,
  file = "../../3_final_figures/LotusSynCom_DA_with_Nodule.rds"
)

ggsave(
  "LotusSynCom_DA_noNodule.pdf",
  plot = main_plot_no_nod,
  width = 21,
  height = 20,
  units = "cm"
)
saveRDS(main_plot_no_nod, file = "LotusSynCom_DA_noNodule.rds")
saveRDS(
  main_plot_no_nod,
  file = "../../3_final_figures/LotusSynCom_DA_noNodule.rds"
)

saveRDS(p_RA, file = "p_RA_Lj_withNodule.rds")
saveRDS(p_RA_no_nod, file = "p_RA_Lj_noNodule.rds")
saveRDS(p_tax_clean, file = "p_tax_clean_Lj.rds")
saveRDS(p_bubble, file = "p_bubble_Lj.rds")
