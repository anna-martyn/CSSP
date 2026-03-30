# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading data
pkg <- c("dplyr", "tidyr", "tibble", "ggplot2", "Maaslin2", "ggtext")
for(pk in pkg){
  library(pk, character.only = T)
}

# Loading data
design <- read.table(
  "../1_data/LotusSC_metadata.txt",
  header = TRUE,
  sep = "\t"
)
asv_table <- read.table(
  "../1_data/LotusSC_ASVtable_nospike.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = ""
)
taxonomy <- read.table(
  "../1_data/LjSC_taxonomy.txt",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)

# Keep only ASVs matched to SynCom
asv_table_matched <- asv_table[grepl("Lj", rownames(asv_table)), ]

# Check if all matched ASVs are present in taxonomy
missing_asvs <- setdiff(rownames(asv_table_matched), taxonomy$ASVid)
if(length(missing_asvs) == 0){
  message("All matched ASVs are present in taxonomy.")
} else {
  warning(length(missing_asvs), " matched ASVs are missing in taxonomy:")
  print(missing_asvs)
}

# Converting ASV reads to relative abundances and adding order level taxonomy
asv_table_norm <- sweep(asv_table_matched, 2, colSums(asv_table_matched), "/")
df <- as.data.frame(asv_table_norm) %>%
  rownames_to_column(var = "ASVid") %>%
  left_join(taxonomy %>% select(ASVid, order), by = "ASVid")

# Long format and adding design
df_long <- df %>%
  pivot_longer(
    cols = -c(ASVid, order),
    names_to = "SampleID",
    values_to = "RA"
  ) %>%
  left_join(design %>% select(SampleID, Compartment, Genotype), by = "SampleID")

df_summary <- df_long %>%
  group_by(SampleID, order, Genotype, Compartment) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop") %>%
  group_by(order, Genotype, Compartment) %>%
  summarise(
    mean_RA = mean(RA, na.rm = TRUE),
    sd = sd(RA, na.rm = T),
    N = length(order),
    .groups = "drop"
  )
df_summary$se <- df_summary$sd/sqrt(df_summary$N)

# Setting compartment order
df_summary <- df_summary %>%
  mutate(
    Compartment = factor(
      Compartment,
      levels = c("Rhizosphere", "Root", "Nodules")
    )
  )

# Lading order colours
colors <- read.table(
  "../../../0_files/Bacterial_order_colors.csv",
  header = TRUE,
  sep = ",",
  comment.char = ""
)

# Setting genotype order with mutant names in italic
df_summary$Genotype <- factor(
  df_summary$Genotype,
  levels = rev(c("WT", "symrk", "ccamk", "nsp1", "nsp2"))
)

genotype_labels <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(size = 6, color = "black"),
  legend.text = element_text(size = 6, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 6),
  legend.position = "right",
  legend.background = element_blank(),
  text = element_text(family = "sans", size = 6, color = "black")
)

# p1 <- ggplot(df_summary, aes(x=Genotype, y=mean_RA, fill=order)) +
#   geom_bar(stat="identity", width=0.5) +
#   facet_wrap(~Compartment, scales="free_x", space = "free_x") +
#   scale_fill_manual(values=colors) +
#   scale_x_discrete(labels = genotype_labels) +
#   scale_y_continuous(expand=c(0,0)) +
#   main_theme +
#   ylab("Mean relative abundance") +
#   labs(fill="Bacterial order") +
#   theme(axis.text.x = element_markdown(size = 6, color="black", angle=30, hjust=1),
#         strip.text.x=element_text(size = 6, face="bold"),
#         legend.position = "bottom",
#         legend.title.position = "top",
#         legend.key.size = unit(0.25, "cm"),
#         legend.key.spacing.y = unit(0, 'cm'),
#         axis.title.x=element_blank()) +
#   guides(fill=guide_legend(nrow=3))+
#   NULL

bar_plot <- ggplot(df_summary, aes(y = Genotype, x = mean_RA, fill = order)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_wrap(
    ~Compartment,
    scales = "free_y",
    space = "free_y",
    nrow = 3,
    strip.position = "right"
  ) +
  scale_fill_manual(values = colors$Color, breaks = colors$Order) +
  scale_x_continuous(expand = c(0, 0)) +
  # scale_y_discrete(labels = genotype_labels, position = "right") +
  scale_y_discrete(labels = genotype_labels) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill = "Bacterial order") +
  theme(
    axis.text.x = element_markdown(size = 6, color = "black"),
    # strip.text.x=element_text(size = 6, face="bold"),
    # strip.text.y.left = element_text(angle = 0, size = 6, face = "bold"),
    strip.text.y = element_text(angle = 0, size = 6, face = "bold"),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.key.size = unit(0.25, "cm"),
    legend.key.spacing.y = unit(0, 'cm'),
    axis.title.x = element_blank(),
    axis.text.y = element_markdown()
  ) +
  guides(fill = guide_legend(nrow = 3)) +
  NULL

# Saving plot
ggsave(
  filename = "1_rds_files/LotusSC_order_RA_stackedbp.pdf",
  plot = bar_plot,
  width = 12,
  height = 6,
  unit = "cm"
)
saveRDS(object = bar_plot, file = "2_figures/LotusSC_order_RA_stackedbp.rds")

# Order level comparisons with Maaslin2 ---------------------------------------
df_summary$Genotype <- factor(
  df_summary$Genotype,
  levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)

# Moving samples and ASV IDs to row names
design <- data.frame(design[,-1], row.names = design$SampleID)
taxonomy <- data.frame(taxonomy[,-1], row.names = taxonomy$ASVid)
orders <- taxonomy[rownames(asv_table_matched), "order"]

# Aggregating ASV table to order level
order_table <- rowsum(asv_table_matched, orders)

# Testing mutant effects for each compartment separately
comps <- c("Rhizosphere", "Root")
summary_list <- list()
for(i in 1:2){
  current_compartment <- comps[i]
  
  # Subsetting data
  samples_sub <- rownames(design)[design$Compartment == current_compartment]
  order_table_sub <- order_table[,samples_sub]
  design_sub <- design[colnames(order_table_sub),]
  
  # Setting factor levels
  design_sub$Genotype <- factor(
    design_sub$Genotype, levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )

  # Differential abundance analysis
  order_test_sub <- Maaslin2(
    input_data = order_table_sub,
    input_metadata = design_sub,
    fixed_effects = "Genotype",
    output = "Maslin_output",
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  
  # Summarising results
  res_sub <- order_test_sub$results
  res_sub <- res_sub[res_sub$qval<0.05, c("feature", "value")]
  colnames(res_sub) <- c("order", "Genotype")
  res_sub$Sig <- "*"
  res_sub$Compartment <- current_compartment

  # Adding DAA results to long form summary table
  df_ord_sub <- df_summary[df_summary$Compartment == current_compartment, ]
  df_ord_sub <- merge(
    df_ord_sub,
    res_sub,
    all.x = TRUE,
    by = c("order", "Compartment", "Genotype")
  )
  df_ord_sub$Sig[is.na(df_ord_sub$Sig)] <- ""
  sub_sig_ord <- unique(df_ord_sub$order[df_ord_sub$Sig == "*"])
  df_ord_sub <- df_ord_sub[df_ord_sub$order %in% sub_sig_ord, ]
  summary_list[[i]] <- df_ord_sub
}

df_ord <- rbind(summary_list[[1]], summary_list[[2]])

# Genotype colours
colors_geno <- c(
  "WT"     = "#A9C289",
  "symrk"  = "#FEDA8B",
  "ccamk"  = "#FDB366",
  "nsp1"   = "#C0E4EF",
  "nsp2"   = "#6EA6CD"
)

# Labels for mutants in italics
genotype_labels_legend <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

df_ord$Genotype <- factor(
  df_ord$Genotype,
  levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)

p_sig <- ggplot(df_ord, aes(x = order, y = mean_RA, fill = Genotype)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    width = 0.8,
    alpha = 0.9
  ) +
  geom_errorbar(
    aes(ymin = mean_RA - 1.96 * se, ymax = mean_RA + 1.96 * se),
    width = 0.3,
    position = position_dodge(width = 0.9)
  ) +
  geom_text(
    data = df_ord,
    aes(
      x = order,
      y = mean_RA + 1.96 * se + 0.015,
      label = Sig,
      fill = Genotype
    ),
    position = position_dodge(width = 0.9),
    inherit.aes = FALSE,
    size = 6
  ) +
  facet_wrap(~Compartment, scales = "free_x") +
  scale_fill_manual(values = colors_geno, labels = genotype_labels_legend) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "", y = "Relative Abundance", fill = "Genotype") +
  main_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.text = element_markdown(),
    strip.text = element_text(face = "bold", size = rel(1)),
    legend.position = "none",
    plot.title = element_text(size = 6)
  )

colnames(df_ord)[1] <- "Order"
df_ord$Host <- "Lotus"

write.csv(x = df_ord, file = "3_tables/df_ord_Lj.csv")

ggsave(
  filename = "2_figures/Lotus_order_RA_sign_orders_asterisks.pdf",
  plot = p_sig,
  width = 12,
  height = 6,
  unit = "cm"
)
saveRDS(
  object = p_sig,
  file = "1_rds_files/Lotus_order_RA_sign_orders_asterisks.rds"
)
