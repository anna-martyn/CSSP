# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
pkg <- c("dplyr", "tidyr", "ggplot2", "tibble", "ggh4x", "scales")
for(pk in pkg){
  library(pk, character.only = T)
}

# Load the Lotus and Hordeum input files.
Lotus_design <- read.table("../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", header=T, sep="\t")
Lotus_asv_table <- read.table( "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t", header = TRUE, row.names = 1, 
  check.names = FALSE, comment.char = "")
Lotus_taxonomy <- read.table("../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv", sep="\t", header=TRUE, fill=TRUE)

Hordeum_design <- read.table("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt", header=T, sep="\t")
Hordeum_asv_table <- read.table(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t", header = TRUE, row.names = 1, 
  check.names = FALSE, comment.char = "", skip = 1)
Hordeum_taxonomy <- read.table("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv",
  sep="\t", header=TRUE, fill=TRUE)

# Clean up the taxonomy file layouts.
rename_tax <- function(tax_table){
  colnames(tax_table)[colnames(tax_table) == "Feature.ID"] <- "ASVid"
  tax_table %>%
    separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
             sep = "; ", fill = "right") %>%
    mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
    replace(is.na(.), "Unknown") %>%
    select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)
}

Lotus_taxonomy <- rename_tax(Lotus_taxonomy)
Hordeum_taxonomy <- rename_tax(Hordeum_taxonomy)

# Re-order the data matrices (only keep samples in design that are found in ASV table, and then filter for ones present in taxonomy file and adjust).
Lotus_design <- Lotus_design %>%
  filter(SampleID %in% colnames(Lotus_asv_table)) 
Lotus_asv_table <- Lotus_asv_table %>% 
  select(all_of(Lotus_design$SampleID)) %>%
  filter(rownames(.) %in% Lotus_taxonomy$ASVid)

Hordeum_design <- Hordeum_design %>%
  filter(SampleID %in% colnames(Hordeum_asv_table)) 
Hordeum_asv_table <- Hordeum_asv_table %>%
  select(all_of(Hordeum_design$SampleID)) %>%
  filter(rownames(.) %in% Hordeum_taxonomy$ASVid)

# Convert ASV reads to relative abundances and save as new dataframe.
asv_to_df <- function(asv_table, taxonomy){
  df <- sweep(asv_table, 2, colSums(asv_table), "/") %>%
    as.data.frame() %>%
    rownames_to_column(var="ASVid") %>%
    left_join(taxonomy %>% select(ASVid, Order), by="ASVid")
  df
}

Lotus_df <- asv_to_df(Lotus_asv_table, Lotus_taxonomy)
Hordeum_df <- asv_to_df(Hordeum_asv_table, Hordeum_taxonomy)

# Reshape to long format and filter for WT and mutant samples in the rhizosphere and root compartment only.
Lotus_df.long <- Lotus_df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(Lotus_design %>% select(SampleID, Soil, Genotype, Compartment), by=c("sampleID"="SampleID")) %>%
  filter(
    Genotype %in% c("WT", "symrk", "ccamk", "nsp1", "nsp2"),
    Compartment %in% c("Rhizosphere", "Root")
  )

Hordeum_df.long <- Hordeum_df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(Hordeum_design %>% select(SampleID, Plant, Soil, Genotype, Compartment), by=c("sampleID"="SampleID")) %>%
  filter(
    Genotype %in% c("WT", "symrk", "ccamk", "nsp1", "nsp2"),
    Compartment %in% c("Rhizosphere", "Root")
  )

# Next we want to choose which bacterial orders to use for representation of the bacterial community structure in a stacked barplot.
# For this, we will first make a summary of the relative abundances by order per sample-soil-compartment group for Lotus and Hordeum.
Lotus_order_summary <- Lotus_df.long %>%
  group_by(Order, sampleID, Soil, Compartment, Genotype) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

Hordeum_order_summary <- Hordeum_df.long %>%
  group_by(Order, sampleID, Soil, Compartment, Genotype) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

## Next, we look at the bacterial orders with highest relative abundances in the samples in order to choose which orders to display.
## We previously did this for the stacked bp showing the differences between Lotus and Hordeum WT across the three soils in both compartments.
## We will do the same here, calculating the top 20 orders in Lotus and Hordeum WT across all soils, and using these orders for displaying the differences
## of those top20 orders in the mutants compared to WT. Also, we filter for the orders with a relative abundance of min. 1% as done previously.

### Step 1: Subset for WT samples.
Lotus_RR <- Lotus_order_summary %>% filter(Genotype == "WT")
Hordeum_RR <- Hordeum_order_summary %>% filter(Genotype == "WT")

### Step 2: Calculate the overall mean RA per order across these samples, identify the top 20 orders by overall mean RA, and keep only those with mean RA ≥ 0.01.
get_top_orders <- function(df){
  df %>%
    group_by(Order) %>%
    summarise(MeanRA = mean(RA, na.rm=T), .groups="drop") %>%
    arrange(desc(MeanRA)) %>%
    slice_head(n=20) %>%
    # filter(MeanRA >= 0.01) %>%
    pull(Order)
}

Lotus_top <- get_top_orders(Lotus_RR)
Hordeum_top <- get_top_orders(Hordeum_RR)
combined_top_orders <- unique(c(Lotus_top, Hordeum_top))
saveRDS(combined_top_orders, "../7_DA_analysis/Orders_to_display.rds")

### Step 3: Combine the orders for Lotus and Hordeum and show how many unique orders would be displayed.
combined_top_orders <- unique(c(Lotus_top, Hordeum_top))
length(combined_top_orders)

# We will work with the selected bacterial orders, go back to the original dataframes and rename all other bacterial orders "Other".
Lotus_df.long <- Lotus_df.long %>%
  mutate(Order = if_else(Order %in% combined_top_orders, Order, "Other"))

Hordeum_df.long <- Hordeum_df.long %>%
  mutate(Order = if_else(Order %in% combined_top_orders, Order, "Other"))

# We then add a new column "Plant" to both dataframes and combine them.
Lotus_df.long <- Lotus_df.long %>% mutate(Plant = "Lotus")
Hordeum_df.long <- Hordeum_df.long %>% mutate(Plant = "Hordeum")

combined_df <- bind_rows(Lotus_df.long, Hordeum_df.long)

# We now summarise the mean RA per plant-genotype-compartment-soil combination.
df.sample_order <- combined_df %>%
  group_by(sampleID, Plant, Compartment, Soil, Genotype, Order) %>%
  summarise(RA=sum(RA), .groups="drop")

df.mean_order <- df.sample_order %>%
  group_by(Plant, Compartment, Soil, Genotype, Order) %>%
  summarise(RA=mean(RA), .groups="drop") %>%
  mutate(Order=factor(Order, levels=c(sort(unique(Order[Order!="Other"])), "Other")),
         Plant=factor(Plant, levels=c("Lotus","Hordeum")),
         Compartment=factor(Compartment, levels=c("Rhizosphere","Root","Nodules")),
         Soil=factor(Soil, levels=c("NPK","PK","UF")),
         Genotype=factor(Genotype, levels=c("WT","symrk","ccamk","nsp1","nsp2")))

# Set the main theme for the stacked barplot.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  panel.border=element_rect(colour="black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(size=8, color="black"),
  legend.text=element_text(size=8, color="black"),
  legend.key=element_blank(),
  axis.title.y=element_text(size=8),
  legend.position="right",
  legend.background=element_blank(),
  text=element_text(family="sans", size=8, color="black")
)

# Set the colours for the bacterial orders.
colors <- read.table("../../../0_files/Bacterial_order_colors.csv", header = T, sep = ",", comment.char = "")

# Make genotype names italic for mutants.
genotype_labels <- c(
  "WT"    = "WT",
  "symrk" = "italic(symrk)",
  "ccamk" = "italic(ccamk)",
  "nsp1"  = "italic(nsp1)",
  "nsp2"  = "italic(nsp2)"
)

# Make the stacked barplot.
p1 <- ggplot(df.mean_order, aes(x=Genotype, y=RA, fill=Order)) +
  geom_bar(stat="identity", width=0.7) +
  scale_fill_manual(values = colors$Color, breaks = colors$Order) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(labels = function(x) parse(text = genotype_labels[x])) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill="Bacterial order") +
  xlab("") +
  guides(fill=guide_legend(nrow=27)) +
  facet_nested(~ Plant + Compartment + Soil, scales="free_x", space="free_x") +
  theme(
    axis.text.x = element_text(size=8, color="black", angle=90, vjust=1, hjust=1),
    strip.text.x = element_text(size=8, face="bold"),
    legend.key.size = unit(0.25, 'cm'),
    legend.margin = margin(l = -8)
  )

p1

# Save the plot.
ggsave("LotusHordeum_Askov_stackedbp_top20_meanRA.pdf", p1, width = 21, height = 6, unit = "cm")
saveRDS(p1, file = "LotusHordeum_Askov_stackedbp_top20_meanRA.rds")
saveRDS(p1, file = "../8_final_figures/LotusHordeum_Askov_stackedbp_top20_meanRA.rds")

# Make a heatmap that displays these data and includes info on significance (differences between WT and individual mutants in each plant-soil-compartment combination.)
## Define input dataframe, filter out Unknown and Other orders, and set factor levels.
df.plot <- df.sample_order %>%
  filter(!(Order %in% c("Other", "Unknown"))) %>%
  mutate(
    Plant = factor(Plant, levels = c("Lotus", "Hordeum")),
    Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
    Genotype = factor(Genotype, levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")),
    Order = factor(Order, levels = sort(unique(Order)))
  )

## Perform significance analysis using linear regression to obtain tests for mutants vs. WT.
mutants <- c("symrk", "ccamk", "nsp1", "nsp2")

Opt <- expand.grid(
  Plant = unique(df.plot$Plant),
  Order = unique(df.plot$Order),
  Compartment = unique(df.plot$Compartment),
  Soil = unique(df.plot$Soil)
)

df_list <- list()
for(i in 1:nrow(Opt)){
  df <- df.plot %>%
  filter(
    Plant == Opt$Plant[i] & 
      Compartment == Opt$Compartment[i] & 
      Soil == Opt$Soil[i] & 
      Order == Opt$Order[i]
  )
  p_vals <- coef(summary(lm(RA~Genotype, data = df)))[-1,"Pr(>|t|)"]
  df_list[[i]] <- tibble(
    Plant = Opt$Plant[i],
    Compartment = Opt$Compartment[i], 
    Soil = Opt$Soil[i],
    Order = Opt$Order[i],
    Genotype = gsub("Genotype", "", names(p_vals)),
    p.value = p_vals
  )
}
df_pvals <- bind_rows(df_list)

# df_pvals <- df.plot %>%
#   group_by(Plant, Compartment, Soil, Order) %>%
#   summarise(
#     p.value = wilcox.test(
#       RA,
#       df.plot$RA[
#         df.plot$Plant == unique(Plant) &
#           df.plot$Compartment == unique(Compartment) &
#           df.plot$Soil == unique(Soil) &
#           df.plot$Order == unique(Order) &
#           df.plot$Genotype == "WT"
#       ],
#       exact = FALSE
#     )$p.value,
#     .groups = "drop"
#   )

## Adjust p-values using Benjamini-Hochberg.
df_pvals <- df_pvals %>%
  group_by(Plant, Compartment, Soil, Order) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(sig = ifelse(p.adj < 0.05, "*", ""))

## Join the significance info with the heatmap data.
df.plot <- df.plot %>%
  left_join(df_pvals %>% select(-p.value), 
            by = c("Plant", "Compartment", "Soil", "Order", "Genotype"))

## Define the heatmap colours and breaks.
breaks <- c(0, 0.005, 0.052, 0.052001, 0.15999, 0.16, 0.34, 0.64)
heat_colors <- c("#1F78B4", "#A6CEE3", "white","#FFFF99",
                 "#FF7F00", "#FB9A99", "#E31A1C", "#902121")
values <- rescale(breaks, to = c(0, 1))

## Make mutant genotype names italic.
genotype_labels_heatmap <- c(
  "WT"    = "WT",
  "symrk" = "italic(symrk)",
  "ccamk" = "italic(ccamk)",
  "nsp1"  = "italic(nsp1)",
  "nsp2"  = "italic(nsp2)"
)

## Make the heatmap plot.
df.plot$Genotype <- factor(df.plot$Genotype, levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2"))
p_heatmap <- ggplot(df.plot, aes(x = Genotype, y = Order, fill = RA)) +
  geom_tile(color = "grey50") +
  geom_text(aes(label = sig), na.rm = TRUE, size = 3) +
  scale_fill_gradientn(colors = heat_colors, values = values,
                       limits = c(0, max(df.plot$RA, na.rm = TRUE)),
                       name = "Relative abundance") +
  scale_y_discrete(limits = rev(levels(df.plot$Order)), position = "right") +
  scale_x_discrete(labels = function(x) parse(text = genotype_labels_heatmap[x])) +
  facet_nested(~ Plant + Compartment + Soil, scales = "free_x", space = "free_x") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(size = 8, angle = 90, vjust = 1, hjust = 1),
    axis.title.y = element_text(size = 8, color = "black"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 8, colour = "black"),
    legend.position = "bottom",
    plot.margin = margin(r = 10, l = 20)
  ) +
  xlab(NULL) +
  ylab("Bacterial order")

p_heatmap

# Save the heatmap.
ggsave("LotusHordeum_Askov_orders_heatmap.pdf", p_heatmap, width = 12, height = 6, unit = "cm")
saveRDS(p_heatmap, file = "LotusHordeum_Askov_orders_heatmap.rds")
saveRDS(p_heatmap, file = "../8_final_figures/LotusHordeum_Askov_orders_heatmap.rds")
