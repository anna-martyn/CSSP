# Seup ------------------------------------------------------------------------
pkg <- c("data.table", "ggplot2")
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Genotype colours
cols <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD",
  "control" = "#cecece"
)

legend_labels <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*",
  "control" = "control"
)

## Loading Lotus data ---------------------------------------------------------
# Lotus feature table filtered for background features
metabolite_data_Lj <- fread(
  "../2_background_removal/1_tables/feature_table_Lotus_filtered.csv"
)

design_Lj <- fread(
  "../1_data/1_Lotus/LotusCSSP_rootex_metadata.txt",
  drop = 4:6
)

# Setting sample names in metadata
design_Lj[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Lj <- design_Lj[ Sample_ID %in% colnames(metabolite_data_Lj)]

# Removing control samples from metadata and feature table
non_control_samples <- design_Lj[Genotype != "control", Sample_ID]
design_Lj <- design_Lj[Sample_ID %in% non_control_samples]
metabolite_data_Lj <- metabolite_data_Lj[,
  c("Feature", non_control_samples),
  with = FALSE
]

## Loading Hordeum data -------------------------------------------------------
metabolite_data_Hv <- fread(
  "../2_background_removal/1_tables/feature_table_Hordeum_filtered.csv"
)
design_Hv <- fread(
  "../1_data/2_Hordeum/HordeumCSSP_rootex_metadata.txt",
  drop = c(2, 4:7)
)

# Setting sample names in metadata
design_Hv[, Sample_ID := paste0("Sample", Sample_ID)]

# Removing samples not in feature table from metadata
design_Hv <- design_Hv[ Sample_ID %in% colnames(metabolite_data_Hv)]

# Adding plant variable to hordeum metadata
design_Hv[,Plant := "Hordeum"]

# Bubble plot ------------------------------------------------------------------
## Lotus -----------------------------------------------------------------------
# Loading results of Genotype effect tests
res_Lj <- fread(
  "../4_genotype_effects/1_tables/Lotus_metabolite_test_results_tobit.csv"
)

# Specifying needed columns
cols_keep <- c(
  "Feature", "symrk", "ccamk", "nsp1", "nsp2",
  "ClassyFire#class", "ClassyFire#class Probability"
)

sum_tab_Lj <- res_Lj[,..cols_keep]
setnames(
  sum_tab_Lj,
  c("ClassyFire#class", "ClassyFire#class Probability"),
  c("Class", "Class_prob")
)

# Setting annotation of features without class as "unknown class"
sum_tab_Lj[is.na(Class) | Class == "", Class := "Unknown class"]

# Setting class of features with class probability less than 60% as 'low probability'
sum_tab_Lj[Class_prob<0.6, Class:="Low probability"]

# Removing columns no longer needed
sum_tab_Lj[,":="(Feature = NULL, Class_prob = NULL)]

# Number of enriched and depleted features
sum_tab_Lj_enr <- sum_tab_Lj[,lapply(.SD, function(x) sum(x == 1) ), Class]
sum_tab_Lj_dep <- sum_tab_Lj[,lapply(.SD, function(x) sum(x == -1) ), Class]

# Long form
sum_tab_Lj_enr <- melt(
  sum_tab_Lj_enr,
  id.vars = "Class",
  variable.name = "Genotype",
  value.name = "Amount_enr"
)
sum_tab_Lj_dep <- melt(
  sum_tab_Lj_dep,
  id.vars = "Class",
  variable.name = "Genotype",
  value.name = "Amount_dep"
)

# Merging enriched and depleted amounts
sum_tab_Lj <- merge(sum_tab_Lj_enr, sum_tab_Lj_dep, by = c("Class", "Genotype"))
sum_tab_Lj[, Host := "Lotus"]

# Bubble plot ------------------------------------------------------------------
## Hordeum ---------------------------------------------------------------------
# Loading results of Genotype effect tests
res_Hv <- fread(
  "../4_genotype_effects/1_tables/Hordeum_metabolite_test_results_tobit.csv"
)

# Specifying needed columns
cols_keep <- c(
  "Feature", "symrk", "ccamk", "nsp1", "nsp2",
  "ClassyFire#class", "ClassyFire#class Probability"
)

sum_tab_Hv <- res_Hv[,..cols_keep]
setnames(
  sum_tab_Hv,
  c("ClassyFire#class", "ClassyFire#class Probability"),
  c("Class", "Class_prob")
)

# Setting annotation of features without class as "unknown class"
sum_tab_Hv[is.na(Class) | Class == "", Class := "Unknown class"]

# Setting class of features with class probability less than 60% as 'low probability'
sum_tab_Hv[Class_prob < 0.6, Class := "Low probability"]

# Removing columns no longer needed
sum_tab_Hv[, ":="(Feature = NULL, Class_prob = NULL)]

# Number of enriched and depleted features
sum_tab_Hv_enr <- sum_tab_Hv[,lapply(.SD, function(x) sum(x == 1) ), Class]
sum_tab_Hv_dep <- sum_tab_Hv[,lapply(.SD, function(x) sum(x == -1) ), Class]

# Long form
sum_tab_Hv_enr <- melt(
  sum_tab_Hv_enr,
  id.vars = "Class",
  variable.name = "Genotype",
  value.name = "Amount_enr"
)
sum_tab_Hv_dep <- melt(
  sum_tab_Hv_dep,
  id.vars = "Class",
  variable.name = "Genotype",
  value.name = "Amount_dep"
)

# Merging enriched and depleted amounts
sum_tab_Hv <- merge(sum_tab_Hv_enr, sum_tab_Hv_dep, by = c("Class", "Genotype"))
sum_tab_Hv[, Host := "Hordeum"]

## Combining tables
sum_tab <- rbind(sum_tab_Lj, sum_tab_Hv)
sum_tab[, Amount := Amount_enr + Amount_dep]

# Total amount of differentially exuded metabolites (DEMs)
total_dem <- sum_tab[,.(N = sum(Amount)), list(Host, Genotype)]

# Total DEMs by class across conditions
total_by_class <- sum_tab[,.(DEMs = max(Amount)), Class]

# Selecting classes to display
classes_to_display <- total_by_class[DEMs > 1, Class]
sum_tab <- sum_tab[Class %in% classes_to_display]

# Setting factor levels
sum_tab[, Host := factor(Host, levels = c("Lotus", "Hordeum"))]

# Removing "unnknown class"
sum_tab <- sum_tab[Class != "Unknown class"]

# Removing rows representing classes with zero DEMs
sum_tab <- sum_tab[Amount != 0]

# Sorting classes by average amount of DEMs across genotypes
avg_sum <- sum_tab[, .(avg = mean(Amount)), Class]
avg_sum <- avg_sum[order(avg)]
sum_tab[,Class:=factor(Class, levels = avg_sum$Class)]

# Creating dummy table used in figure to make room for the total DEMs and 
# adding to sum_tab
dummy <- data.table(
  Class = "Total DEMs",
  Genotype = rep(c("symrk", "ccamk", "nsp1", "nsp2"), 2),
  Amount_enr = NA,
  Amount_dep = NA,
  Host = rep(c("Lotus", "Hordeum"), each = 2),
  Amount = NA
)
sum_tab <- rbind(sum_tab, dummy)

# Adding the proportion of enriched DEMs by class and genotype
sum_tab[, enriched_proportion := Amount_enr / Amount]

# Table containing text to be printed on figure
text_dt <- data.table(
  Class = "Total DEMs",
  Host = rep(c("Lotus", "Hordeum"), each = 4),
  Genotype = rep(c("symrk", "ccamk", "nsp1", "nsp2"), 2)
)
text_dt <- merge(text_dt, total_dem, by = c("Host", "Genotype"))

# Setting factor levels
sum_tab[,Host:=factor(Host, levels = c("Lotus", "Hordeum"))]

# Plot
bubble_plot <- ggplot(
  data = sum_tab,
  mapping = aes(
    x = Genotype,
    y = Class,
    size = Amount,
    fill = enriched_proportion
  )
) +
  geom_label(
    data = text_dt,
    mapping = aes(x = Genotype, y = Class, label = N),
    size = 6 / .pt,
    fill = NA,
    linewidth = NA
  ) +
  geom_point(shape = 21) +
  facet_wrap(~ factor(Host, levels = c("Lotus", "Hordeum"))) +
  scale_fill_gradient2(
    midpoint = 0.5,
    low = "darkblue",
    mid = "white",
    high = "#902121",
    name = "Proportion of enriched DEMs"
  ) +
  scale_size_continuous(breaks = c(1, 10, 50, 100), range = c(1, 4)) +
  labs(size = "# of DEMs") +
  ggtitle("") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    strip.background = element_rect(colour = NA),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 6,
      colour = "black"
    ),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6, colour = "black"),
    strip.text = element_text(colour = 'black', size = 6, face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.justification = c(1.5, 0),
    legend.text = element_text(size = 6, margin = margin(t = 2, unit = "pt")),
    legend.title = element_text(size = 6, margin = margin(r = 5)),
    legend.key = element_rect(fill = NA),
    legend.key.spacing.x = unit(5, "pt"),
    plot.title = element_text(size = 6, hjust = 0.5),
    legend.margin = margin(t = 2, r = 8, l = 5, unit = "pt")
  ) +
  scale_x_discrete(
    labels = c(
      "symrk" = expression(italic("symrk")),
      "ccamk" = expression(italic("ccamk")),
      "nsp1" = expression(italic("nsp1")),
      "nsp2" = expression(italic("nsp2"))
    )
  ) +
  guides(fill = guide_colourbar(barheight = 0.5, label.position = "bottom")) +
  NULL

# Siving plot
saveRDS(object = bubble_plot, file = "1_rds_files/bubble_plot.rds")
ggsave(
  filename = "2_figures/bubble_plot.pdf",
  plot = bubble_plot,
  height = 18,
  width = 9,
  units = "cm"
)

