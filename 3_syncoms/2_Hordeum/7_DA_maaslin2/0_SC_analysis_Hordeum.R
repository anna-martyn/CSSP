## Original script by Ib Thorsgaard Jensen (Aalborg University), modified and extended by Anna Martyn (amartyn@mpipz.mpg.de)

# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Loading packages
pkg <- c(
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "ggh4x",
  "ggpubr",
  "cowplot",
  "patchwork",
  "Maaslin2"
)
for(pk in pkg) library(pk, character.only = TRUE)

# Setting working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading file with structural zeros function
source("Structural_zeros.R")

# Loading data
design <- read.table(
  "../1_data/HordeumSC_metadata.txt",
  header = TRUE,
  sep = "\t"
)
asv_table <- read.table(
  "../1_data/HordeumSC_ASVtable.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1
)
taxonomy <- read.table(
  "../1_data/CerealSC_taxonomy_May23.txt",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)

# Keep only matched ASVs in ASV table
asv_table_matched <- asv_table[grepl("_", rownames(asv_table)), , drop = FALSE]

# Library sizes for each sample (needed for Structural_zeros3)
design$library_size <- colSums(asv_table_matched[, design$SampleID])

# Splitting data by compartment
samples_rhizo <- design$SampleID[design$Compartment == "Rhizosphere"]
samples_root  <- design$SampleID[design$Compartment == "Root"]

asv_table_rhizo <- asv_table_matched[, samples_rhizo, drop = FALSE]
asv_table_root <- asv_table_matched[, samples_root, drop = FALSE]

meta_rhizo <- design %>%
  filter(SampleID %in% samples_rhizo) %>%
  column_to_rownames("SampleID")

meta_root <- design %>%
  filter(SampleID %in% samples_root) %>%
  column_to_rownames("SampleID")

# Genotype factor levels
meta_rhizo$Genotype <- factor(
  meta_rhizo$Genotype,
  levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)
meta_root$Genotype <- factor(
  meta_root$Genotype,
  levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
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
# Rhizosphere
da_sz_rhizo <- Structural_zeros3(
  asv_table_rhizo,
  meta_rhizo,
  group = "Genotype",
  ref = "WT",
  min_reads = 20,
  min_present_reps = 2
)

# Root
da_sz_root <- Structural_zeros3(
  asv_table_root,
  meta_root,
  group = "Genotype",
  ref = "WT",
  min_reads = 20,
  min_present_reps = 2
)

# Differential abundance analysis ----------------------------------------------
# Rhizosphere
da_maa_rhizo <- Maaslin2(
  input_data = asv_table_rhizo,
  input_metadata = meta_rhizo,
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
struc_zero_wt <- da_sz_rhizo$struc_zero_table[all_da_rhizo,1]
struc_zero_mutants <- da_sz_rhizo$struc_zero_table[all_da_rhizo,-1]
da_sz_sign_rhizo <- struc_zero_wt - struc_zero_mutants

for(g in c("symrk", "ccamk", "nsp1", "nsp2")){
  rows <- da_sz_rhizo$struc_zero_DA[[g]]
  cols <- paste(c("DA", "Lfc"), g, sep = "_")
  results_rhizo[rows, cols[1]] <- da_sz_sign_rhizo[rows, g]
  results_rhizo[rows, cols[2]] <- da_sz_sign_rhizo[rows, g] * Inf
}

# Root
da_maa_root <- Maaslin2(
  input_data = asv_table_root,
  input_metadata = meta_root,
  output = "Maaslin2_root",
  fixed_effects = "Genotype",
  min_prevalence = 0.1,
  plot_heatmap = F,
  plot_scatter = F
)
res <- da_maa_root$results

for(g in c("symrk", "ccamk", "nsp1", "nsp2")){
  res_g <- res[res$value == g, ]
  cols <- paste(c("Lfc", "DA"), g, sep = "_")
  rows <- res_g$feature
  results_root[rows, cols[1]] <- res_g$coef
  results_root[rows, cols[2]] <- (res_g$qval < 0.05) * sign(res_g$coef)
}

# Updating results with structural zeros
all_da_root <- Reduce("union", da_sz_root$struc_zero_DA)
struc_zero_wt <- da_sz_root$struc_zero_table[all_da_root, 1]
struc_zero_mutants <- da_sz_root$struc_zero_table[all_da_root,-1]
da_sz_sign_root <- struc_zero_wt - struc_zero_mutants

for(g in c("symrk","ccamk","nsp1","nsp2")){
  rows <- da_sz_root$struc_zero_DA[[g]]
  cols <- paste(c("DA", "Lfc"), g, sep = "_")
  results_root[rows, cols[1]]  <- da_sz_sign_root[rows, g]
  results_root[rows, cols[2]] <- da_sz_sign_root[rows, g] * Inf
}

# Merging results with taxonomy and relative abundance info
RA_rhizo <- t(t(asv_table_rhizo)/colSums(asv_table_rhizo))
RA_root  <- t(t(asv_table_root)/colSums(asv_table_root))

# Rhizosphere
results_rhizo <- merge(
  x = results_rhizo,
  y = taxonomy,
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

# Root
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
write.csv(x = results_rhizo, file = "3_tables/DA_SynCom_Hordeum_rhizo.csv")
write.csv(x = results_root, file = "3_tables/DA_SynCom_Hordeum_root.csv")

# Plot highlighting DA ASVs ----------------------------------------------------

# Mean RA for WT in Rhizosphere and Root
asv_table_RA <- sweep(asv_table_matched, 2, colSums(asv_table_matched), "/")
asv_RA_long <- asv_table_RA %>%
  rownames_to_column("ASVid") %>%
  pivot_longer(cols = -ASVid, names_to = "SampleID", values_to = "RA") %>%
  left_join(design %>% select(SampleID, Compartment, Genotype), by = "SampleID")

# Keep only WT samples
asv_RA_wt <- asv_RA_long %>%
  filter(Genotype == "WT") %>%
  group_by(ASVid, Compartment) %>%
  summarise(mean_RA = mean(RA, na.rm = TRUE), .groups = "drop")

# Adding taxonomy
asv_RA_wt <- asv_RA_wt %>%
  left_join(taxonomy %>% select(ASVid, Order), by = "ASVid") %>%
  mutate(Order = ifelse(is.na(Order), "Unknown", Order))

# Ordering ASVs by taxonomic order
asv_order_levels <- asv_RA_wt %>%
  distinct(ASVid, Order) %>%
  arrange(Order) %>%
  pull(ASVid)

asv_RA_wt$ASVid <- factor(asv_RA_wt$ASVid, levels = asv_order_levels)

# Compartment factor levels
asv_RA_wt$Compartment <- factor(
  asv_RA_wt$Compartment,
  levels = c("Rhizosphere", "Root")
)

## Taxonomic color bar ---------------------------------------------------------
colors <- read.table(
  "../../../0_files/Bacterial_order_colors.csv",
  header = TRUE,
  sep = ",",
  comment.char = ""
)

tax_bar <- asv_RA_wt %>%
  distinct(ASVid, Order) %>%
  filter(ASVid %in% isolate_keep)

p_tax <- ggplot(tax_bar, aes(x = ASVid, y = 1, fill = Order)) +
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
    plot.margin = margin(0.25, 0, 0.25, 0, unit = "lines")
  ) +
  # guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))+
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
  labs(y = "Differential abundance\nin mutants") +
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
    axis.title.y = element_blank(),
    strip.placement = "outside",
    legend.position = "bottom",
    axis.text.y = element_blank(),
    legend.text = element_text(color = "black", size = 6),
    legend.title = element_text(color = "black", size = 6),
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "lines"),
    panel.spacing = unit(0.4, "lines")
  ) +
  NULL

## Mean RA bar plots -----------------------------------------------------------
### All compartments -----------------------------------------------------------

# Changing names for appropriate line breaks in plot
asv_RA_wt$Compartment <- as.character(asv_RA_wt$Compartment)
asv_RA_wt$Compartment[
  asv_RA_wt$Compartment == "Rhizosphere"
] <- "Rhizo-\nsphere"
asv_RA_wt$Compartment <- factor(
  asv_RA_wt$Compartment,
  levels = c("Rhizo-\nsphere", "Root")
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
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ggtitle("Hordeum") +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(color = "black", size = 6, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "lines"),
    panel.spacing = unit(0.4, "lines")
  ) +
  force_panelsizes(cols = c(1, 1, 1), rows = c(1, 1, 0.7)) +
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
  ) +
  NULL

## Combining plots -------------------------------------------------------------
# Remove individual legends from plots
p_tax_clean <- p_tax + theme(legend.position = "none")

# Combine plots vertically
main_plot <- p_RA /
  p_tax_clean /
  (p_bubble + theme(legend.position = "none")) +
  plot_layout(heights = c(0.35, 0.05, 0.6))

main_plot <- main_plot +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.spacing = unit(0, "cm"))

lgd_bubble <- ggpubr::get_legend(p_bubble, position = "bottom")
lgd_tax <- ggpubr::get_legend(p_tax, position = "bottom")
lgd <- plot_grid(lgd_bubble, lgd_tax, ncol = 1)

# Saving plot
ggsave(
  filename = "2_figures/HordeumSynCom_DA_noNodule.pdf",
  plot = main_plot,
  width = 21,
  height = 20,
  units = "cm"
)
saveRDS(object = main_plot, file = "1_rds_files/HordeumSynCom_DA_noNodule.rds")

saveRDS(object = p_RA, file = "1_rds_files/p_RA_Hv_no_nodule.rds")
saveRDS(object = p_tax_clean, file = "1_rds_files/p_tax_clean_Hv.rds")
saveRDS(object = p_bubble, file = "1_rds_files/p_bubble_Hv.rds")
saveRDS(object = lgd, file = "1_rds_files/legend.rds")

