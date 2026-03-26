# Genotype-specific differential abundance analysis using MaAsLin2 and
# structural zero detection, followed by visualisation.

# Setup ------------------------------------------------------------------------
# Clean up
options(warn = -1)
rm(list = ls())

# Loading packages
pkg <- c("data.table", "ggplot2", "Maaslin2", "ggh4x", "ggtext")

for(pk in pkg){
  library(pk, character.only = T)
}

# Setting directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load file with function for structural zero analysis
source("../0_files/Structural_zeros.R")

# Lading data
asv_table <- fread(
  "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv"
)
colnames(asv_table)[1] <- "ASVid"

design <- fread(
  "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt",
  drop = c(5, 7, 8)
)

taxonomy <- fread(
  "../../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv"
)

# Cleaning up taxonomy
setnames(taxonomy, "Feature ID", "ASVid")
taxa_levels <- c(
    "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
  )
taxonomy[,c(taxa_levels):= tstrsplit(Taxon, "; ", fill = "u__Unknown")]
taxonomy[,c(taxa_levels):= lapply(.SD, substr, 4, 1000), .SDcols = taxa_levels]
taxonomy[Kingdom == "ssigned", Kingdom:= "Unassigned"]
taxonomy[,Taxon:=NULL]
setcolorder(taxonomy, c("ASVid", taxa_levels, "Confidence"))

# Merging ASV table and metadata
asv_table_with_metadata <- transpose(
  asv_table,
  make.names = "ASVid",
  keep.names = "SampleID"
)
asv_table_with_metadata <- merge(
  design,
  asv_table_with_metadata,
  by = "SampleID"
)

# Vector of mutants
mutants <- c("symrk", "ccamk", "nsp1", "nsp2")

# Setting reference levels
asv_table_with_metadata[,
  Genotype := factor(
    Genotype,
    levels = c("WT", mutants)
  )
]

opt <- as.matrix(
  expand.grid(
    Compartment = c("Root", "Rhizosphere"),
    Soil = c("NPK", "PK", "UF")
  )
)

# Differential abundance analysis with MaAslin2 --------------------------------
# Subsetting data for each compartment-soil combination
asv_tables_subsets <- apply(
  opt, 1,
  function(x){
    # Subsetting data for compartment-soil combination
    data_subset <- asv_table_with_metadata[Compartment == x[1] & Soil == x[2]]
    
    # Extracting metadata and dropping unused factor levels
    col_design <- c("Plant", "Genotype", "Compartment", "Soil")
    design <- data.frame(
      data_subset[, ..col_design],
      row.names = data_subset$SampleID
    )
    design$Genotype <- droplevels(design$Genotype)
    
    # Extracting and transposing ASV table for use in MaAsLin2
    asv_table <- t(data_subset[,-(1:5)])
    colnames(asv_table) <- data_subset$SampleID
    
    # Returning ASV and metadata as list
    out <- list(asv_table = asv_table, design = design)
    return(out)
  }
)

# Running MaAsLin2 and structural zero analysis for each subset
daa_results <- lapply(asv_tables_subsets, function(X){
  # Creating output directory if it doesn't already exist
  if(!dir.exists("Maaslin2 files")){
    dir.create("Maaslin2 files")
  }
  
  # MaAsLin2
  da_maaslin <- Maaslin2(
    input_data = X$asv_table,
    input_metadata = X$design,
    output = "Maaslin2 files",
    max_significance = 0.05,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    fixed_effects = "Genotype",
    min_prevalence = 0.1,
    reference = "WT"
  )
  res_dt <- data.table(da_maaslin$results)
  res_dt[, feature := gsub("X", "", feature)]
  res_dt <- res_dt[order(feature)]

  # Structural zeros
  X$design$library_size <- colSums(X$asv_table)

  da_struc_zero <- Structural_zeros3(
    X$asv_table,
    X$design,
    group = "Genotype",
    ref = "WT",
    min_reads = 20,
    min_present_reps = 3
  )

  # Constructing binary matrix to indicate which ASVs are DA
  res_mat <- matrix(0, nrow = nrow(X$asv_table), ncol = 4)
  rownames(res_mat) <- sort(rownames(X$asv_table))
  colnames(res_mat) <- mutants

  for(i in 1:4){
    current_mutant <- mutants[i]
    da_wt_vs_mutant <- res_dt[value == current_mutant & qval < 0.05]$feature
    res_mat[da_wt_vs_mutant, current_mutant] <- sign(
      res_dt[value == current_mutant & qval < 0.05]$coef
    )
  }
  
  # Summarising ASVs DA by structural zero
  da_asv_struc_zero_all <- unlist(da_struc_zero$struc_zero_DA)
  da_asv_amount_by_genotype <- unlist(lapply(da_struc_zero$struc_zero_DA, length))
  names(da_asv_struc_zero_all) <- rep(
    mutants,
    da_asv_amount_by_genotype
  )

  struc_zero_res_dt <- data.table(
    ASV = da_asv_struc_zero_all,
    Genotype = names(da_asv_struc_zero_all),
    DA_by_struc_zero = TRUE
  )

  # Filling in structural zeros results
  for(i in 1:4){
    current_mutant <- mutants[i]
    struc_zero_wt_vs_mutant <- setdiff(
      da_struc_zero$struc_zero_DA[current_mutant][[1]],
      res_dt[value == current_mutant & qval < 0.05]$feature
    )
    da_struc_zero_wt_vs_mutant <- da_struc_zero$struc_zero_table[
      struc_zero_wt_vs_mutant,
      c("WT", current_mutant)
    ]
    struc_zero_signs <- apply(
      da_struc_zero_wt_vs_mutant[, c(current_mutant, "WT")],
      1,
      diff
    )
    res_mat[names(struc_zero_signs), current_mutant] <- struc_zero_signs
  }

  # Mean RA of DA ASVs in WT samples
  wt_samples <- rownames(X$design)[X$design$Genotype == "WT"]
  lib_size <- colSums(X$asv_table[, wt_samples])
  mean_RAs <- apply(
    res_mat,
    2,
    function(y) {
      asvs <- rownames(res_mat)[y != 0]
      asv_subset <- X$asv_table[asvs, wt_samples]
      mean(rowSums(t(asv_subset) / lib_size))
    }
  )

  mean_RA_res <- data.table(
    Compartment = X$design$Compartment[1],
    Soil = X$design$Soil[1],
    Genotype = names(mean_RAs),
    mean_RAs
  )

  # Summarising and combining results
  res_dt[, ":="(
    DA_by_sig = qval < 0.05,
    metadata = NULL,
    stderr = NULL,
    pval = NULL,
    N = NULL,
    N.not.zero = NULL,
    name = NULL,
    Compartment = X$design$Compartment[1],
    Soil = X$design$Soil[1]
  )]
  colnames(res_dt)[1:3] <- c("ASV", "Genotype", "LogFC")
  res_dt <- merge(
    res_dt,
    struc_zero_res_dt,
    by = c("ASV", "Genotype"),
    all.x = TRUE
  )
  res_dt[is.na(DA_by_struc_zero), DA_by_struc_zero := FALSE]
  res_dt[, qval := NULL]
  res_dt <- res_dt[DA_by_sig | DA_by_struc_zero]

  # Number of DA ASVs per condition
  amount <- res_dt[,
    .(DA_ASVs = sum(DA_by_sig | DA_by_struc_zero)),
    list(Genotype)
  ]
  amount[, ":="(
    Compartment = X$design$Compartment[1],
    Soil = X$design$Soil[1],
    TOtal_ASVs = nrow(da_maaslin$results)/4
  )]
  setcolorder(
    amount,
    c("Genotype", "Compartment", "Soil", "TOtal_ASVs", "DA_ASVs")
  )

  out <- list(
    res_mat = res_mat,
    res_dt = res_dt,
    amount = amount,
    mean_RA_res = mean_RA_res
  )

  return(out)
})

# Exporting summary tables
da_asv_amount <- rbindlist(
  lapply(daa_results, function(x) x$amount)
)
da_asv_amount[, ":="(
  Genotype = factor(Genotype, levels = mutants),
  Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
  Soil = factor(Soil, levels = c("NPK", "PK", "UF"))
)]
da_asv_amount <- da_asv_amount[order(Soil, Compartment, Genotype)]
fwrite(da_asv_amount, "Lotus_DA_ASVs_overview.csv")

da_asv_res <- rbindlist(
  lapply(daa_results, function(x) x$res_dt)
)
setcolorder(da_asv_res, c("ASV", "Soil", "Compartment", "Genotype"))
da_asv_res[,
  Genotype := factor(
    Genotype,
    levels = c("WT", mutants)
  )
]
da_asv_res <- da_asv_res[order(Soil, Compartment, Genotype, ASV)]
fwrite(da_asv_res, "Lotus_DAanalysis_results.csv")

mean_RA_res <- rbindlist(
  lapply(daa_results, function(x) x$mean_RA_res)
)

# Preparing results for heatmap ------------------------------------------------
# Combining DA matrices across all conditions
res_mat_full <- t(Reduce("cbind", lapply(daa_results, function(x) x$res_mat)))
res_mat_full <- data.table(Genotype = rownames(res_mat_full), res_mat_full)

# Long format for plot
heatmap_data <- lapply(
  daa_results,
  function(X){
    dt <- data.table(ASVid = rownames(X$res_mat), X$res_mat)
    dtt <- melt(dt, id.vars = 1, variable.name = "Genotype", value.name = "DAA")
    dtt[, ":="(Compartment = X$res_dt$Compartment[1], Soil = X$res_dt$Soil[1])]
    return(dtt)
  }
)
heatmap_data <- rbindlist(heatmap_data)

# Converting numeric DA indicators to categorical variable
heatmap_data[,
  DAA := fcase(
    DAA == 0  , "NS"       ,
    DAA == 1  , "Enriched" ,
    DAA == -1 , "Depleted"
  )
]

# Keeping only DA ASVs in heatmap
da_asv <- unique(heatmap_data[DAA != "NS"]$ASVid)
htmp_hiabn <- heatmap_data[ASVid %in% da_asv]

# Taxonomy annotation for heatmap ----------------------------------------------
# Subsetting taxonomy to include only DA ASVs
tax_bar <- taxonomy[ASVid %in% unique(htmp_hiabn$ASVid)]

# Load order colours
colors_orders <- fread("../../../../0_files/Bacterial_order_colors.csv")

# Loading orders to display. These are top 20 mean RA orders in either Lotus or
# Hordeum WT across compartments and soils (as used for supplementary figures)
combined_top_orders <- readRDS("../Orders_to_display.rds")

# Assigning orders without assigned colours to 'Other'
tax_bar[!(Order %in% combined_top_orders), Order := "Other"]
tax_bar[, Order := factor(Order, levels = colors_orders$Order)]
tax_bar <- tax_bar[order(Order)]

# Matching ASV ordering between taxonomy and heatmap
htmp_hiabn <- htmp_hiabn[,ASVid:=factor(ASVid, levels = tax_bar$ASVid)]
htmp_hiabn <- htmp_hiabn[order(ASVid)]
tax_bar$ASVid <- factor(tax_bar$ASVid, levels = tax_bar$ASVid)

# Taxonomy bar
p_tax <- ggplot(tax_bar, aes(x = ASVid, y = 1, fill = Order)) +
  geom_tile() +
  scale_fill_manual(
    values = colors_orders$Color,
    breaks = colors_orders$Order,
    drop = FALSE
  ) +
  theme_void() +
  labs(fill = "Bacterial order") +
  theme(
    legend.position = "none",
    legend.text = element_text(color = "black", size = 6),
    legend.title = element_text(color = "black", size = 6, face = "bold"),
    legend.key.size = unit(0.25, 'cm'),
    legend.key.spacing.y = unit(0, 'cm'),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "lines")
  ) +
  NULL

# Heatmap for differential abundance --------------------------------------
# Merging DA counts for y-axis annotation
htmp_hiabn <- merge(
  htmp_hiabn, da_asv_amount, by = c("Compartment", "Soil", "Genotype")
)

# Formatting y-axis labels to include number of DA ASVs
htmp_hiabn[,
  y_axis := paste0("*", Genotype, "*", " (", DA_ASVs, ")")
]

# Genotype ordering
ordering <- c(
  unique(htmp_hiabn[grepl("nsp2", y_axis)]$y_axis),
  unique(htmp_hiabn[grepl("nsp1", y_axis)]$y_axis),
  unique(htmp_hiabn[grepl("ccamk", y_axis)]$y_axis),
  unique(htmp_hiabn[grepl("symrk", y_axis)]$y_axis)
)

htmp_hiabn[, y_axis := factor(y_axis, levels = ordering)]

# Heatmap
heatmap <- ggplot(
  data = htmp_hiabn,
  mapping = aes(x = ASVid, y = y_axis, fill = DAA)
) +
  geom_tile() +
  facet_wrap2(
    vars(Compartment, Soil),
    strip = strip_nested(),
    ncol = 1,
    strip.position = "left",
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c("darkblue", "#902121", "white"),
    breaks = c("Depleted", "Enriched", "NS")
  ) +
  labs(y = "Lotus") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.2, 'lines'),
    strip.background = element_rect(colour = NA),
    strip.placement = "outside",
    axis.title.y = element_text(
      size = 6,
      family = "Helvetica",
      colour = "black",
      face = "bold"
    ),
    axis.text.y = element_markdown(
      size = 6,
      family = "Helvetica",
      colour = "black"
    ),
    legend.text = element_text(size = 6, family = "Helvetica"),
    legend.title = element_text(size = 6, family = "Helvetica"),
    strip.text = element_text(size = 6, family = "Helvetica", face = "bold"),
    plot.margin = margin(t = 0, r = 0, b = 0.5, l = 0, unit = "lines")
  ) +
  NULL

# Barplot of cumulative relative abundance of DA ASVs in WT --------------------
mean_RA_res[,
  Genotype := factor(Genotype, levels = rev(mutants))
]
bar_plot <- ggplot(
  data = mean_RA_res,
  mapping = aes(x = mean_RAs, y = Genotype)
) +
  geom_bar(stat = "identity") +
  facet_wrap2(
    vars(Compartment, Soil),
    strip = strip_nested(),
    ncol = 1,
    strip.position = "left"
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.1, 0.2)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    panel.spacing = unit(0.2, 'lines'),
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 6,
      family = "Helvetica",
      colour = "black"
    ),
    legend.text = element_text(size = 6, family = "Helvetica"),
    legend.title = element_text(size = 6, family = "Helvetica"),
    strip.text = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0.5, l = 0, unit = "lines")
  ) +
  NULL

# Axis title panel
p_axis_title <- ggplot() +
  theme_void() +
  geom_text(aes(x = -10, y = 0, label = "RA in WT"),
            fontface = "bold", size = 6/.pt) +
  xlab(NULL) + ylab(NULL)

# Saving plots -----------------------------------------------------------------
saveRDS(heatmap, "LotusCSSP_Askov_DA_heatmap.rds")
saveRDS(bar_plot, "LotusCSSP_Askov_DA_barplot.rds")
saveRDS(p_tax, "LotusCSSP_Askov_DA_taxonomy.rds")
saveRDS(p_axis_title, "LotusCSSP_Askov_DA_axis_title.rds")
