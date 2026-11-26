# Lotus SynCom experiment - bacterial load calculation

# Background info --------------------------------------------------------------
# During library preparation, we added 10^4 copies of synthetic spike-in DNA to
# PCR1 of each Lotus SynCom sample. As the same number of spike-in copies was
# added to every sample, the ratio of bacterial to spike-in reads reflects the
# amount of bacterial 16S rRNA template per PCR. We use this ratio to estimate
# and compare the total and commensal bacterial load of wild-type and mutant samples.

# No spike-in was used during Hordeum library preparation, hence we do not have
# data for Hordeum samples.

# Spike-in ASVs were identified in the ASV table clean-up script. In the SynCom
# dataset, only the spike-in ASV itself (85fa8bb918a926d97659d9b64ca6fedd) was
# found and there is no spike-in variant.

# As in all other SynCom analyses, only ASVs assigned to SynCom members are
# counted as bacterial reads. Before removing ASVs not matching any SynCom member,
# we however inspect their share and test whether it differs between wild-type and mutants.

# The script consists of a check of non-SynCom ASVs, followed by two parts:
# Part 1: Spike-in inspection and total bacterial load
# Part 2: Nodule ASVs, symbiont/nodule-ASV definition and commensal bacterial load

# Set working directory --------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load required packages -------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)

# Settings ---------------------------------------------------------------------

# Input and output files -------------------------------------------------------
feature_file  <- "../1_data/LotusSC_ASVtable.tsv"
metadata_file <- "../1_data/LotusSC_metadata.txt"
taxonomy_file <- "../1_data/LjSC_taxonomy.txt"
out_dir       <- "Results"
dir.create(out_dir, showWarnings = FALSE)

# Spike-in info ----------------------------------------------------------------
spike_ids    <- c("85fa8bb918a926d97659d9b64ca6fedd")
spike_copies <- 1e4

# Symbiont/nodule-ASV definition ----------------------------------------------------------
# To compare commensal levels, we need to define which ASVs are symbionts/nodule-ASVs.
# We use two definitions:
# 1) ASVs of the order Rhizobiales with a mean relative abundance (RA) of min. 1%
# in wild-type nodules, and
# 2) all ASVs with min. 1% mean RA in wild-type nodules (including non-rhizobial
# nodule co-colonisers). All remaining ASVs are considered commensals.
nodule_threshold <- 0.01
symbiont_order   <- "Rhizobiales"

# To test whether the commensal load results depend on the chosen threshold, we
# repeat the comparison using thresholds of 0.5% and 5%.
sensitivity_thresholds <- c(0.005, 0.01, 0.05)

# Factor levels and colours ----------------------------------------------------
genotype_order    <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
compartment_order <- c("Rhizosphere", "Root", "Nodules")
colors <- c("WT" = "#A9C289", "symrk" = "#FEDA8B", "ccamk" = "#FDB366",
            "nsp1" = "#C0E4EF", "nsp2" = "#6EA6CD")

# Functions --------------------------------------------------------------------

# Mutant vs. wild-type comparison ----------------------------------------------
# Wilcoxon rank-sum test of each mutant vs. wild-type, separately for each compartment.
# P-values are adjusted for multiple testing using Benjamini-Hochberg.
# We use: fold change = median of mutant / median of wild-type.
compare_vs_wt <- function(df) {
  df %>%
    group_by(compartment) %>%
    group_modify(function(d, key) {
      bind_rows(lapply(setdiff(levels(droplevels(d$genotype)), "WT"), function(m) {
        x <- d$value[d$genotype == "WT"]
        y <- d$value[d$genotype == m]
        if (length(x) < 2 || length(y) < 2) return(NULL)
        tibble(mutant = m, n_wt = length(x), n_mut = length(y),
               fold_change = median(y) / median(x),
               p = wilcox.test(y, x, exact = FALSE)$p.value)
      }))
    }) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p, method = "BH"))
}

# Significance labels from adjusted p-values
stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                               p < 0.05 ~ "*", TRUE ~ "ns")

# Print statistics
print_stats <- function(stats) {
  stats %>%
    mutate(fold_change = round(fold_change, 2), p_adj = signif(p_adj, 2),
           sig = stars(p_adj)) %>%
    select(compartment, mutant, n_wt, n_mut, fold_change, p_adj, sig) %>%
    arrange(compartment, mutant) %>%
    as.data.frame() %>% print(row.names = FALSE)
}

# Define scientific labels format
sci_labels <- function(x) {
  parse(text = gsub("e\\+?", " %*% 10^", scales::scientific_format()(x)))
}

# Make plots -------------------------------------------------------------------
# Boxplot per genotype, faceted by compartment using a log10 y-axis
base_plot <- function(d, ylab) {
  ggplot(d, aes(genotype, value, fill = genotype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.2) +
    facet_wrap(~ compartment, scales = "free_y", nrow = 1) +
    scale_y_log10(breaks = scales::breaks_log(n = 5), labels = sci_labels) +
    scale_fill_manual(values = colors) +
    labs(x = NULL, y = ylab) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# Add significance labels
add_stars <- function(p, d, stats) {
  ann <- d %>%
    group_by(compartment) %>%
    summarise(y = max(value) * 2, .groups = "drop") %>%
    right_join(stats, by = "compartment") %>%
    mutate(genotype = factor(mutant, levels = genotype_order),
           label = stars(p_adj))
  p + geom_text(data = ann, aes(genotype, y, label = label),
                inherit.aes = FALSE, size = 3.2)
}

# Save plots
save_plot <- function(p, name, width = 4, height = 3) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height)
  saveRDS(p, file.path(out_dir, paste0(name, ".rds")))
}

# Load data --------------------------------------------------------------------
feat <- read_tsv(feature_file, skip = 1, show_col_types = FALSE)
colnames(feat)[1] <- "ASV"
stopifnot(all(spike_ids %in% feat$ASV))

# Taxonomy of the SynCom members
tax <- read_tsv(taxonomy_file, show_col_types = FALSE)
colnames(tax)[1] <- "ASV"
tax <- tax %>% mutate(Taxon = paste(order, genus, sep = "; "))
symbiont_order_ids <- tax$ASV[tax$order == symbiont_order]

meta <- read_tsv(metadata_file, show_col_types = FALSE) %>%
  rename(sample_id = SampleID, genotype = Genotype, compartment = Compartment) %>%
  mutate(genotype    = factor(genotype, levels = genotype_order),
         compartment = factor(compartment, levels = compartment_order))
stopifnot(!any(is.na(meta$genotype)))

# Check whether all samples in metadata are present in ASV table.
missing <- setdiff(meta$sample_id, colnames(feat))
if (length(missing)) print(paste("Not in ASV table:", paste(missing, collapse = ", ")))

# Convert ASV table to long format and only keep samples in metadata.
counts <- feat %>%
  pivot_longer(-ASV, names_to = "sample_id", values_to = "count") %>%
  filter(sample_id %in% meta$sample_id)

# Non-SynCom ASVs --------------------------------------------------------------
# Before removing ASVs not matching any SynCom member, we check how many reads
# they account for, and whether their share differs between wild-type and mutants.

# Share of non-SynCom reads per sample (spike-in reads not counted)
non_syncom <- counts %>%
  filter(!ASV %in% spike_ids) %>%
  group_by(sample_id) %>%
  summarise(non_syncom_pct = 100 * sum(count[!ASV %in% tax$ASV]) / sum(count),
            .groups = "drop") %>%
  left_join(meta, by = "sample_id")
write_tsv(non_syncom, file.path(out_dir, "Non_SynCom_share_per_sample.tsv"))

# Overview per compartment
non_syncom %>%
  group_by(compartment) %>%
  summarise(pct_median = round(median(non_syncom_pct), 2),
            pct_min = round(min(non_syncom_pct), 2),
            pct_max = round(max(non_syncom_pct), 2), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# Mutant vs. wild-type (nodules excluded, as only in wild-type)
stats_non_syncom <- compare_vs_wt(non_syncom %>%
                                    filter(compartment != "Nodules") %>%
                                    mutate(value = non_syncom_pct))
write_tsv(stats_non_syncom, file.path(out_dir, "Non_SynCom_share_stats.tsv"))
print_stats(stats_non_syncom)

# Plot share of non-SynCom reads.
d0 <- non_syncom %>% filter(compartment != "Nodules") %>% mutate(value = non_syncom_pct)
ann0 <- d0 %>%
  group_by(compartment) %>%
  summarise(y = max(value) * 1.1, .groups = "drop") %>%
  right_join(stats_non_syncom, by = "compartment") %>%
  mutate(genotype = factor(mutant, levels = genotype_order), label = stars(p_adj))
p0 <- ggplot(d0, aes(genotype, value, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2) +
  geom_text(data = ann0, aes(genotype, y, label = label),
            inherit.aes = FALSE, size = 3.2) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Non-SynCom reads (%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p0, "0_non_syncom_share")
print(p0)

# Remove non-SynCom ASVs and keep SynCom members and spike-in.
counts <- counts %>% filter(ASV %in% c(tax$ASV, spike_ids))

# Part 1: Spike-in inspection and total bacterial load -------------------------

# Spike-in reads per sample ----------------------------------------------------
# For each sample, we sum the reads of the spike-in and of all bacterial (SynCom)
# ASVs, and calculate the percentage of spike-in reads.
per_sample <- counts %>%
  group_by(sample_id) %>%
  summarise(spike_reads     = sum(count[ASV %in% spike_ids]),
            bacterial_reads = sum(count[!ASV %in% spike_ids]), .groups = "drop") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(spike_pct = round(100 * spike_reads / (spike_reads + bacterial_reads), 1))

# Make an overview of spike-in reads and percentage for each compartment.
per_sample %>%
  group_by(compartment) %>%
  summarise(reads_median = median(spike_reads), reads_min = min(spike_reads),
            reads_max = max(spike_reads), pct_median = median(spike_pct),
            pct_min = min(spike_pct), pct_max = max(spike_pct), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# Inspect samples sorted by spike-in percentage.
per_sample %>%
  arrange(desc(spike_pct)) %>%
  select(sample_id, compartment, genotype, spike_reads, bacterial_reads, spike_pct) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Inspect samples sorted by number of bacterial reads.
per_sample %>%
  arrange(bacterial_reads) %>%
  select(sample_id, compartment, genotype, spike_reads, bacterial_reads, spike_pct) %>%
  as.data.frame() %>% print(row.names = FALSE)

# No sample was dominated by spike-in reads and therefore no sample had to be excluded.
excluded_samples <- c()
per_sample <- per_sample %>% mutate(included = !sample_id %in% excluded_samples)

# Make first plot showing spike-in reads per sample
pA <- base_plot(per_sample %>%
                  filter(compartment != "Nodules") %>%
                  mutate(value = spike_reads),
                "Spike-in reads per sample")
save_plot(pA, "1_spike_reads")
print(pA)

# Total bacterial load ---------------------------------------------------------
per_sample <- per_sample %>%
  mutate(total_load = bacterial_reads / spike_reads * spike_copies)

# Mutant vs. wild-type
# Note: We do not include nodule samples here as only wild-type had nodules.
cmp <- per_sample %>% filter(included, compartment != "Nodules")

stats_total <- compare_vs_wt(cmp %>% mutate(value = total_load))
write_tsv(stats_total, file.path(out_dir, "Total_load_stats.tsv"))
print_stats(stats_total)

# Make second plot showing total bacterial load.
dB <- cmp %>% mutate(value = total_load)
pB <- add_stars(base_plot(dB, "Total bacterial load"),
                dB, stats_total)
save_plot(pB, "2_total_bacterial_load")
print(pB)

# Part 2: Nodule ASVs, symbionts and commensal bacterial load ------------------

# Relative abundance of bacterial ASVs -----------------------------------------
# Remove spike-in reads and add sample info.
bact <- counts %>%
  filter(!ASV %in% spike_ids) %>%
  group_by(sample_id) %>% mutate(ra = count / sum(count)) %>% ungroup() %>%
  left_join(meta %>% select(sample_id, compartment, genotype), by = "sample_id")

# Nodule ASVs ------------------------------------------------------------------
# Calculate mean RA of each ASV in wild-type nodules.
nod_ra <- bact %>%
  filter(compartment == "Nodules", genotype == "WT") %>%
  group_by(ASV) %>%
  summarise(mean_ra_wt_nodules = mean(ra),
            n_nodule_samples   = sum(ra > 0), .groups = "drop")

# Define all nodule ASVs with min. 1% mean RA in wild-type nodules,
# including ALL bacterial orders.
nod_ids <- nod_ra$ASV[nod_ra$mean_ra_wt_nodules >= nodule_threshold]

# Inspect the mean RA of these ASVs in all compartments for each genotype.
nodule_asv_check <- bact %>%
  filter(ASV %in% nod_ids) %>%
  group_by(ASV, compartment, genotype) %>%
  summarise(mean_ra = round(mean(ra), 3), .groups = "drop") %>%
  left_join(tax %>% select(ASV, Taxon), by = "ASV") %>%
  mutate(rhizobiales = ASV %in% symbiont_order_ids) %>%
  pivot_wider(names_from = genotype, values_from = mean_ra) %>%
  group_by(ASV) %>%
  mutate(wt_nodule_ra = WT[compartment == "Nodules"]) %>%
  ungroup() %>%
  arrange(desc(wt_nodule_ra), ASV, compartment) %>%
  select(ASV, Taxon, rhizobiales, compartment, all_of(genotype_order))
write_tsv(nodule_asv_check, file.path(out_dir, "Nodule_asvs_RA_by_compartment.tsv"))

nodule_asv_check %>% as.data.frame() %>% print(row.names = FALSE)

# Symbiont/nodule-ASV definitions ---------------------------------------------------------
# We use two definitions of symbionts/nodule-ASVs (= non-commensal ASVs):

# 1) "rhizobiales" -> All Rhizobiales ASVs with min. mean RA in wild-type nodules.
# Only includes rhizobial symbionts and treats other ASVs as commensal.

# 2) "all_nodule" -> All ASVs with min. mean RA in wild-type nodules included.
# Here we include all nodule colonisers.

symbiont_definitions <- c("rhizobiales", "all_nodule")

symbiont_asvs <- function(th, definition) {
  s <- nod_ra$ASV[nod_ra$mean_ra_wt_nodules >= th]
  if (definition == "rhizobiales") s <- intersect(s, symbiont_order_ids)
  s
}

# Calculate the symbiont/nodule-ASV and commensal share per sample.
read_shares <- function(th, definition) {
  bact %>%
    group_by(sample_id) %>%
    summarise(symbiont_share  = sum(ra[ASV %in% symbiont_asvs(th, definition)]),
              commensal_share = 1 - symbiont_share, .groups = "drop")
}

# Save the list of symbiont/nodule ASVs for both definitions.
bind_rows(lapply(symbiont_definitions, function(def) {
  nod_ra %>%
    filter(ASV %in% symbiont_asvs(nodule_threshold, def)) %>%
    mutate(definition = def)
})) %>%
  left_join(tax %>% select(ASV, Taxon), by = "ASV") %>%
  arrange(definition, desc(mean_ra_wt_nodules)) %>%
  write_tsv(file.path(out_dir, "Symbiont_asvs.tsv"))

# Mean symbiont share per compartment and genotype, for both definitions.
# This shows how much commensal load can differ from total load at all.
symbiont_share <- bind_rows(lapply(symbiont_definitions, function(def) {
  read_shares(nodule_threshold, def) %>% mutate(definition = def)
})) %>%
  filter(!sample_id %in% excluded_samples) %>%
  left_join(meta, by = "sample_id") %>%
  group_by(definition, compartment, genotype) %>%
  summarise(mean_symbiont_share = round(mean(symbiont_share), 3), .groups = "drop")
write_tsv(symbiont_share, file.path(out_dir, "Symbiont_share_by_group.tsv"))

symbiont_share %>%
  filter(compartment == "Root") %>%
  pivot_wider(names_from = genotype, values_from = mean_symbiont_share) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Commensal load ---------------------------------------------------------------
# We use: commensal load = total load * share of commensal (non-symbiont/nodule) 
# reads, calculated for both symbiont/nodule-ASV definitions.
per_sample <- per_sample %>%
  left_join(read_shares(nodule_threshold, "rhizobiales") %>%
              select(sample_id, commensal_share_rhizobiales = commensal_share),
            by = "sample_id") %>%
  left_join(read_shares(nodule_threshold, "all_nodule") %>%
              select(sample_id, commensal_share_all_nodule = commensal_share),
            by = "sample_id") %>%
  mutate(commensal_load_rhizobiales = total_load * commensal_share_rhizobiales,
         commensal_load_all_nodule  = total_load * commensal_share_all_nodule)
write_tsv(per_sample, file.path(out_dir, "Load_per_sample.tsv"))

# Mutant vs. wild-type (for both definitions)
cmp <- per_sample %>% filter(included, compartment != "Nodules")

stats_commensal_rhizobiales <- compare_vs_wt(cmp %>% mutate(value = commensal_load_rhizobiales))
write_tsv(stats_commensal_rhizobiales, file.path(out_dir, "Commensal_load_stats_rhizobiales.tsv"))
print_stats(stats_commensal_rhizobiales)

stats_commensal_all_nodule <- compare_vs_wt(cmp %>% mutate(value = commensal_load_all_nodule))
write_tsv(stats_commensal_all_nodule, file.path(out_dir, "Commensal_load_stats_all_nodule.tsv"))
print_stats(stats_commensal_all_nodule)

# Perform a sensitivity analysis
# Does the commensal load comparison change if we change the nodule RA threshold
sens <- bind_rows(lapply(symbiont_definitions, function(def) {
  bind_rows(lapply(sensitivity_thresholds, function(th) {
    cmp %>%
      left_join(read_shares(th, def) %>% select(sample_id, cs = commensal_share),
                by = "sample_id") %>%
      mutate(value = total_load * cs) %>%
      compare_vs_wt() %>%
      mutate(definition = def, threshold = th)
  }))
}))
write_tsv(sens, file.path(out_dir, "Threshold_sensitivity.tsv"))

sens %>%
  group_by(definition, threshold) %>%
  summarise(n_tests = n(), n_significant = sum(p_adj < 0.05),
            n_higher_in_mutant = sum(p_adj < 0.05 & fold_change > 1),
            n_lower_in_mutant  = sum(p_adj < 0.05 & fold_change < 1),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# Make third plot showing commensal load (here symbionts/nodule-ASVs = Rhizobiales nodule ASVs).
dC1 <- cmp %>% mutate(value = commensal_load_rhizobiales)
pC1 <- add_stars(base_plot(dC1, "Commensal load\n(Rhizobiales nodule ASVs removed)"),
                 dC1, stats_commensal_rhizobiales)
save_plot(pC1, "3_commensal_load_rhizobiales")
print(pC1)

# Make fourth plot showing commensal load (here symbionts/nodule-ASVs = all nodule ASVs).
dC2 <- cmp %>% mutate(value = commensal_load_all_nodule)
pC2 <- add_stars(base_plot(dC2, "Commensal load\n(all nodule ASVs removed)"),
                 dC2, stats_commensal_all_nodule)
save_plot(pC2, "4_commensal_load_all_nodule_ASVs")
print(pC2)

# Combine plots into supplementary figure --------------------------------------
# Define the desired text size.
small_text <- theme(
  text = element_text(size = 6),
  axis.text = element_text(size = 6),
  strip.text = element_text(size = 6)
)

text_6pt <- function(p) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomText"))
      p$layers[[i]]$aes_params$size <- 6 / .pt
  }
  p + small_text
}

# Load plots.
spike_syncom     <- readRDS(file.path(out_dir, "1_spike_reads.rds"))
total_syncom     <- readRDS(file.path(out_dir, "2_total_bacterial_load.rds"))
comm_rhiz_syncom <- readRDS(file.path(out_dir, "3_commensal_load_rhizobiales.rds"))
comm_all_syncom  <- readRDS(file.path(out_dir, "4_commensal_load_all_nodule_ASVs.rds"))

# Combine plots.
final_plot_syncom <- plot_grid(
  text_6pt(spike_syncom),
  text_6pt(total_syncom),
  text_6pt(comm_rhiz_syncom),
  text_6pt(comm_all_syncom),
  ncol = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 8,
  align = "hv",
  axis = "tblr"
)

# Save supplementary figure.
ggsave(filename = file.path("Suppl_Figure_spikein_SynCom.pdf"),
       plot = final_plot_syncom, width = 18, height = 12, unit = "cm")
