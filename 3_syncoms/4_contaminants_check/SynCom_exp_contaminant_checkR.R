# Lotus and Hordeum SynCom experiments - non-SynCom ASVs (contaminants)

# Background info --------------------------------------------------------------
# In the SynCom experiments, ASVs were assigned to SynCom members using the 16S
# rRNA V5-V7 sequences of the inoculated isolates (VSEARCH, 99% identity). ASVs
# not matching any SynCom member (e.g. seed-borne or contaminating bacteria) were
# removed prior to downstream analysis.

# Here, we check for each host separately:
# 1) how many reads these non-SynCom ASVs account for in rhizosphere and root,
#    and whether this share differs between wild-type and CSSP mutants;
# 2) whether wild-type and mutants contain the same non-SynCom ASVs;
# 3) for Hordeum, which samples carry a dominant seed-borne Enterobacter ASV, and
#    whether this affects the beta diversity (PCoA) and differential abundance
#    results of the SynCom members.

# In Lotus, spike-in DNA was added during library preparation; spike-in reads
# are removed before calculating the share of non-SynCom reads. No spike-in was
# used for Hordeum.

# Set working directory --------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load required packages -------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(vegan)
library(readxl)
library(cowplot)

# Settings ---------------------------------------------------------------------

# Input files per host ---------------------------------------------------------
cfg_lotus <- list(
  host          = "Lotus",
  feature_file  = "../1_Lotus/1_data/LotusSC_ASVtable.tsv",
  metadata_file = "../1_Lotus/1_data/LotusSC_metadata.txt",
  taxonomy_file = "../1_Lotus/1_data/LjSC_taxonomy.txt",
  spike_ids     = c("85fa8bb918a926d97659d9b64ca6fedd"))

cfg_hordeum <- list(
  host          = "Hordeum",
  feature_file  = "../2_Hordeum/1_data/HordeumSC_ASVtable.tsv",
  metadata_file = "../2_Hordeum/1_data/HordeumSC_metadata.txt",
  taxonomy_file = "../2_Hordeum/1_data/CerealSC_taxonomy_May23.txt",
  spike_ids     = c())

# Differential abundance results of the Hordeum SynCom experiment (MaAsLin2;
# one sheet per compartment)
da_file <- "../2_Hordeum/7_DA_maaslin2/3_tables/DA_SynCom_Hordeum.xlsx"

# Enterobacter ASV detected in part of the Hordeum wild-type plants, and the
# share of reads above which a sample is considered to carry it
entero_id        <- "931ceae7b121c3a084311d252fdda9a5"
entero_threshold <- 5

out_dir <- "Results"
dir.create(out_dir, showWarnings = FALSE)

# Samples with fewer reads are excluded (same criterion as for alpha diversity)
min_reads <- 1000

# Factor levels and colours ----------------------------------------------------
genotype_order    <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
compartment_order <- c("Rhizosphere", "Root")
colors <- c("WT" = "#A9C289", "symrk" = "#FEDA8B", "ccamk" = "#FDB366",
            "nsp1" = "#C0E4EF", "nsp2" = "#6EA6CD")

# Functions --------------------------------------------------------------------

# Load ASV table, metadata and SynCom taxonomy of one host ----------------------
# Returns the counts in long format (rhizosphere and root only, spike-in removed)
# with a column indicating whether an ASV is a SynCom member.
load_host <- function(cfg) {
  feat <- read_tsv(cfg$feature_file, skip = 1, show_col_types = FALSE)
  colnames(feat)[1] <- "ASV"
  
  # All ASVs listed in the taxonomy file are SynCom members
  tax <- read_tsv(cfg$taxonomy_file, show_col_types = FALSE)
  syncom_ids <- tax[[1]]
  
  meta <- read_tsv(cfg$metadata_file, show_col_types = FALSE) %>%
    dplyr::rename(sample_id = SampleID, genotype = Genotype, compartment = Compartment) %>%
    filter(compartment %in% compartment_order) %>%
    mutate(genotype    = factor(genotype, levels = genotype_order),
           compartment = factor(compartment, levels = compartment_order))
  stopifnot(!any(is.na(meta$genotype)))
  
  feat %>%
    pivot_longer(-ASV, names_to = "sample_id", values_to = "count") %>%
    filter(sample_id %in% meta$sample_id, !ASV %in% cfg$spike_ids) %>%
    mutate(syncom = ASV %in% syncom_ids, host = cfg$host) %>%
    left_join(meta %>% dplyr::select(sample_id, genotype, compartment), by = "sample_id")
}

# Mutant vs. wild-type comparison ----------------------------------------------
# Wilcoxon rank-sum test of each mutant vs. wild-type, separately for each
# compartment. P-values are adjusted for multiple testing using Benjamini-Hochberg.
# We use: fold change = median of mutant / median of WT.
compare_vs_wt <- function(df) {
  df %>%
    group_by(compartment) %>%
    group_modify(function(d, key) {
      bind_rows(lapply(setdiff(levels(droplevels(d$genotype)), "WT"), function(m) {
        x <- d$value[d$genotype == "WT"]
        y <- d$value[d$genotype == m]
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


# Load data --------------------------------------------------------------------
counts <- bind_rows(load_host(cfg_lotus), load_host(cfg_hordeum)) %>%
  mutate(host = factor(host, levels = c("Lotus", "Hordeum")))


# Part 1: Share of non-SynCom reads --------------------------------------------

# Percentage of reads per sample that belong to non-SynCom ASVs
non_syncom <- counts %>%
  group_by(host, sample_id, genotype, compartment) %>%
  summarise(total_reads    = sum(count),
            non_syncom_pct = 100 * sum(count[!syncom]) / sum(count),
            .groups = "drop")

# Inspect samples with low read numbers, then exclude them
non_syncom %>% filter(total_reads < min_reads) %>%
  as.data.frame() %>% print(row.names = FALSE)
low_samples <- non_syncom$sample_id[non_syncom$total_reads < min_reads]
non_syncom  <- non_syncom %>% filter(!sample_id %in% low_samples)
counts      <- counts %>% filter(!sample_id %in% low_samples)
write_tsv(non_syncom, file.path(out_dir, "Non_SynCom_share_per_sample.tsv"))

# Overview per host and compartment
non_syncom %>%
  group_by(host, compartment) %>%
  summarise(pct_median = round(median(non_syncom_pct), 1),
            pct_min = round(min(non_syncom_pct), 1),
            pct_max = round(max(non_syncom_pct), 1), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# Mutant vs. wild-type, separately for each host (BH correction within each host)
stats_share <- non_syncom %>%
  mutate(value = non_syncom_pct) %>%
  group_by(host) %>%
  group_modify(function(d, key) compare_vs_wt(d)) %>%
  ungroup()
write_tsv(stats_share, file.path(out_dir, "Non_SynCom_share_stats.tsv"))

stats_share %>%
  mutate(fold_change = round(fold_change, 2), p_adj = signif(p_adj, 2),
         sig = stars(p_adj)) %>%
  dplyr::select(host, compartment, mutant, n_wt, n_mut, fold_change, p_adj, sig) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Plot share of non-SynCom reads per host and compartment
ann <- non_syncom %>%
  group_by(host, compartment) %>%
  summarise(y = max(non_syncom_pct) * 1.1, .groups = "drop") %>%
  right_join(stats_share, by = c("host", "compartment")) %>%
  mutate(genotype = factor(mutant, levels = genotype_order), label = stars(p_adj))

p_share <- ggplot(non_syncom, aes(genotype, non_syncom_pct, fill = genotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2) +
  geom_text(data = ann, aes(genotype, y, label = label),
            inherit.aes = FALSE, size = 3.2) +
  facet_grid(compartment ~ host, scales = "free_y") +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Non-SynCom reads (%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "Non_SynCom_share.pdf"), p_share, width = 5, height = 4.5)
print(p_share)


# Part 2: Do wild-type and mutants contain the same non-SynCom ASVs? -----------

# Non-SynCom reads only
contam <- counts %>% filter(!syncom, count > 0)

# 2a) Shared non-SynCom ASVs ---------------------------------------------------
# For each host, compartment and mutant, we calculate:
#   - the share of non-SynCom reads in the mutant that belong to ASVs also
#     detected in wild-type (in at least one sample), and
#   - the share of non-SynCom reads in wild-type that belong to ASVs also
#     detected in the mutant.
# Reads (rather than numbers of ASVs) are used, as many non-SynCom ASVs are rare
# and only detected sporadically in single samples.
reads_shared <- function(d, from, to) {
  asv_to <- unique(d$ASV[d$genotype == to])
  x <- d[d$genotype == from, ]
  round(100 * sum(x$count[x$ASV %in% asv_to]) / sum(x$count), 1)
}

shared <- contam %>%
  group_by(host, compartment) %>%
  group_modify(function(d, key) {
    bind_rows(lapply(setdiff(genotype_order, "WT"), function(m) {
      tibble(mutant = m,
             mut_reads_in_WT_ASVs_pct = reads_shared(d, m, "WT"),
             WT_reads_in_mut_ASVs_pct = reads_shared(d, "WT", m))
    }))
  }) %>%
  ungroup()
write_tsv(shared, file.path(out_dir, "Non_SynCom_shared_reads.tsv"))
shared %>% as.data.frame() %>% print(row.names = FALSE)

# 2b) Most abundant non-SynCom ASVs --------------------------------------------
# Mean share (% of all reads) of the five most abundant non-SynCom ASVs per host,
# per compartment and genotype. This shows whether single ASVs dominate in
# specific genotypes (e.g. 931ceae7b121c3a084311d252fdda9a5 in Hordeum wild-type,
# which matches Enterobacter (BLAST)).
top_asvs <- contam %>%
  group_by(host, ASV) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  group_by(host) %>%
  slice_max(total, n = 5, with_ties = FALSE) %>%
  ungroup()

counts %>%
  group_by(sample_id) %>%
  mutate(pct = 100 * count / sum(count)) %>%
  ungroup() %>%
  semi_join(top_asvs, by = c("host", "ASV")) %>%
  group_by(host, compartment, ASV = substr(ASV, 1, 8), genotype) %>%
  summarise(mean_pct = round(mean(pct), 2), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = mean_pct) %>%
  arrange(host, compartment, desc(WT)) %>%
  as.data.frame() %>% print(row.names = FALSE)


# Part 3: Enterobacter in Hordeum ----------------------------------------------

# 3a) Share of Enterobacter per sample -----------------------------------------
# Share of Enterobacter reads (% of all reads) per Hordeum sample. Samples above
# the threshold are considered to carry Enterobacter.
entero <- counts %>%
  filter(host == "Hordeum") %>%
  group_by(sample_id, genotype, compartment) %>%
  summarise(entero_pct = 100 * sum(count[ASV == entero_id]) / sum(count),
            .groups = "drop") %>%
  mutate(entero = entero_pct > entero_threshold)
write_tsv(entero, file.path(out_dir, "Hordeum_Enterobacter_per_sample.tsv"))

entero %>%
  filter(entero_pct > 0.5) %>%
  mutate(entero_pct = round(entero_pct, 1)) %>%
  arrange(compartment, genotype, sample_id) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Wild-type samples with and without Enterobacter
wt_entero <- entero$sample_id[entero$genotype == "WT" & entero$entero]
wt_clean  <- entero$sample_id[entero$genotype == "WT" & !entero$entero]

# 3b) PCoA of SynCom members with Enterobacter-carrying samples marked ---------
# As in Fig. 6C: relative abundance of SynCom ASVs only, Bray-Curtis, PCoA per
# compartment.
pcoa <- bind_rows(lapply(compartment_order, function(comp) {
  m <- counts %>%
    filter(host == "Hordeum", syncom, compartment == comp) %>%
    dplyr::select(sample_id, ASV, count) %>%
    pivot_wider(names_from = ASV, values_from = count, values_fill = 0)
  ra  <- as.matrix(m[, -1]) / rowSums(m[, -1])
  res <- cmdscale(vegdist(ra, method = "bray"), k = 2, eig = TRUE)
  tibble(sample_id = m$sample_id, PCo1 = res$points[, 1], PCo2 = res$points[, 2],
         compartment = comp)
})) %>%
  left_join(entero %>% dplyr::select(sample_id, genotype, entero), by = "sample_id") %>%
  mutate(compartment = factor(compartment, levels = compartment_order))

# Position of wild-type samples
pcoa %>% filter(genotype == "WT") %>%
  mutate(across(c(PCo1, PCo2), ~ round(.x, 2))) %>%
  arrange(compartment, PCo1) %>%
  as.data.frame() %>% print(row.names = FALSE)

p_pcoa <- ggplot(pcoa, aes(PCo1, PCo2, colour = genotype, shape = entero)) +
  geom_point(size = 2) +
  facet_wrap(~ compartment, scales = "free") +
  scale_colour_manual(values = colors) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("FALSE" = "no", "TRUE" = "yes"),
                     name = "Enterobacter") +
  labs(title = "Hordeum") +
  theme_bw()
print(p_pcoa)

# 3c) Differential abundance compared with wild-type plants without Enterobacter
# For each isolate with significant differential abundance (DA) in a mutant, we
# check whether the difference also holds when the mutant is compared with
# wild-type plants without Enterobacter only (same direction and at least
# 1.4-fold, i.e. |log2 fold change| > 0.5). With only two such plants, this is a
# descriptive plausibility check, not a statistical test.
ra_hv <- counts %>%
  filter(host == "Hordeum", syncom) %>%
  group_by(sample_id) %>%
  mutate(ra = count / sum(count)) %>%
  ungroup()

da_check <- bind_rows(lapply(compartment_order, function(comp) {
  da <- read_excel(da_file, sheet = comp, na = "NA")
  bind_rows(lapply(setdiff(genotype_order, "WT"), function(m) {
    da %>%
      dplyr::select(ASV = ASVid, Genus, DA = all_of(paste0("DA_", m))) %>%
      mutate(DA = as.numeric(DA)) %>%
      filter(!is.na(DA), DA != 0) %>%
      mutate(compartment = comp, mutant = m)
  }))
})) %>%
  rowwise() %>%
  mutate(
    ra_mut      = mean(ra_hv$ra[ra_hv$ASV == ASV & ra_hv$compartment == compartment &
                                  ra_hv$genotype == mutant]),
    ra_wt_clean = mean(ra_hv$ra[ra_hv$ASV == ASV & ra_hv$compartment == compartment &
                                  ra_hv$sample_id %in% wt_clean]),
    lfc_clean   = log2((ra_mut + 1e-5) / (ra_wt_clean + 1e-5))) %>%
  ungroup() %>%
  mutate(consistent = sign(lfc_clean) == sign(DA) & abs(lfc_clean) > 0.5)
write_tsv(da_check, file.path(out_dir, "Hordeum_DA_check_WT_without_Enterobacter.tsv"))

# Summary: number of DA isolates per mutant and how many are consistent
da_check %>%
  group_by(compartment, mutant) %>%
  summarise(n_DA = n(), n_consistent = sum(consistent), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

# Isolates whose difference does not hold
da_check %>% filter(!consistent) %>%
  dplyr::select(compartment, mutant, ASV, Genus, DA, lfc_clean) %>%
  mutate(lfc_clean = round(lfc_clean, 2)) %>%
  as.data.frame() %>% print(row.names = FALSE)

# Plot: number of DA isolates per mutant, coloured by whether the difference
# also holds compared with wild-type plants without Enterobacter
p_da <- da_check %>%
  mutate(mutant      = factor(mutant, levels = genotype_order),
         compartment = factor(compartment, levels = compartment_order),
         consistent  = factor(ifelse(consistent, "yes", "no"), levels = c("no", "yes"))) %>%
  ggplot(aes(mutant, fill = consistent)) +
  geom_bar() +
  facet_wrap(~ compartment) +
  scale_fill_manual(values = c("yes" = "grey40", "no" = "grey85"),
                    name = "Consistent vs. wild-type\nwithout Enterobacter") +
  labs(x = NULL, y = "Differentially abundant isolates", title = "Hordeum") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_da)


# Supplementary figure ---------------------------------------------------------
p_share_si <- p_share + facet_wrap(~ host + compartment, scales = "free_y",
                                   nrow = 2, dir = "v")

# Same text size (6 pt) for all panels, including plot titles and legends
small_text <- theme(text = element_text(size = 6),
                    axis.text = element_text(size = 6),
                    strip.text = element_text(size = 6),
                    plot.title = element_text(size = 6),
                    legend.text = element_text(size = 6),
                    legend.title = element_text(size = 6))
p_share_si$layers[[3]]$aes_params$size <- 6 / .pt

fig_si <- plot_grid(p_share_si + small_text,
                    plot_grid(p_pcoa + small_text, p_da + small_text,
                              ncol = 2, labels = c("B", "C"), label_size = 8),
                    ncol = 1, labels = c("A", ""), label_size = 8,
                    rel_heights = c(1.2, 1))
ggsave(file.path(out_dir, "Suppl_Figure_contaminants.pdf"), fig_si,
       width = 18, height = 20, units = "cm")
