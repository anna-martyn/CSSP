# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages.
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggh4x)
library(FSA)# Clean up.
options(warn=-1)
library(multcompView)
library(scales)

# Load Lotus and Hordeum files.# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages.
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggh4x)
library(FSA)# Clean up.
options(warn=-1)
library(multcompView)
library(scales)


# Load Lotus and Hordeum files.
Lotus_design <- read.table(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", header=T, sep="\t"
)
Lotus_asv_table <- read.table(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t", header = TRUE, row.names = 1, 
  check.names = FALSE, comment.char = ""
)
Lotus_taxonomy <- read.table(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv",
  sep="\t", header=TRUE, fill=TRUE
)

Hordeum_design <- read.table(
  "../../1_data/2_Barley/HordeumCSSP_AskovSoils_metadata.txt", header=T, sep="\t"
)
Hordeum_asv_table <- read.table(
  "../../1_data/2_Barley/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t", header = TRUE, row.names = 1, 
  check.names = FALSE, comment.char = "", skip = 1
)
Hordeum_taxonomy <- read.table(
  "../../1_data/2_Barley/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv",
  sep="\t", header=TRUE, fill=TRUE
)

# Rename columns Feature.ID to ASVid.
colnames(Lotus_taxonomy)[colnames(Lotus_taxonomy) == "Feature.ID"] <- "ASVid"
colnames(Hordeum_taxonomy)[colnames(Hordeum_taxonomy) == "Feature.ID"] <- "ASVid"

# Clean-up layout of taxonomy files.
Lotus_taxonomy <- Lotus_taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right") %>%
  mutate(across(Kingdom:Species, ~sub("^.{3}", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

Hordeum_taxonomy <- Hordeum_taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right") %>%
  mutate(across(Kingdom:Species, ~sub("^.{3}", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Re-order data matrices.
Lotus_design <- Lotus_design %>%
  filter(SampleID %in% colnames(Lotus_asv_table)) ## Subset design to keep samples present in ASV table.

Lotus_asv_table <- Lotus_asv_table %>% ## Subset and reorder the ASV table to match the filtered design file and only keep ASVs present in the taxonomy file.
  select(all_of(Lotus_design$SampleID)) %>%
  filter(rownames(.) %in% Lotus_taxonomy$ASVid)

Hordeum_design <- Hordeum_design %>%
  filter(SampleID %in% colnames(Hordeum_asv_table)) 

Hordeum_asv_table <- Hordeum_asv_table %>%
  select(all_of(Hordeum_design$SampleID)) %>%
  filter(rownames(.) %in% Hordeum_taxonomy$ASVid)

# Convert ASV reads to relative abundances and save as new dataframe.
Lotus_asv_table <- sweep(Lotus_asv_table, 2, colSums(Lotus_asv_table), "/")
Lotus_df <- as.data.frame(Lotus_asv_table) %>%
  rownames_to_column(var="ASVid") %>%
  left_join(Lotus_taxonomy %>% select(ASVid, Order), by="ASVid")

Hordeum_asv_table <- sweep(Hordeum_asv_table, 2, colSums(Hordeum_asv_table), "/")
Hordeum_df <- as.data.frame(Hordeum_asv_table) %>%
  rownames_to_column(var="ASVid") %>%
  left_join(Hordeum_taxonomy %>% select(ASVid, Order), by="ASVid")

# Reshape to long format and filter for WT samples.
Lotus_df.long <- Lotus_df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(Lotus_design %>% select(SampleID, Soil, Genotype, Compartment), by=c("sampleID"="SampleID")) %>%
  filter(Genotype=="WT")

Hordeum_df.long <- Hordeum_df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(Hordeum_design %>% select(SampleID, Plant, Soil, Genotype, Compartment), by=c("sampleID"="SampleID")) %>%
  filter(Genotype=="WT")

# Next we want to choose which bacterial orders to use for representation of the bacterial community structure in a stacked barplot.

## For this, we will first make a summary of the relative abundances by order per sample - soil - compartment group for Lotus and Hordeum.
Lotus_order_summary <- Lotus_df.long %>%
  group_by(Order, sampleID, Soil, Compartment) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

Hordeum_order_summary <- Hordeum_df.long %>%
  group_by(Order, sampleID, Soil, Compartment) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

## Next, we look at the bacterial orders with highest relative abundances in the samples in order to choose which orders to display.

## There are different ways to do this:

### 1) One can look at the top20 bacterial orders with highest mean relative abundance in the different Soil - Compartment combinations for Lotus and Hordeum separately.
##### With this option one would however end up with a large number of different bacterial orders.
##### (Even when filtering for orders with e.g. mean relative abundance of min. 1%).

### 2) One can go by the top20 bacterial orders in either the rhizosphere OR the root of Lotus or Hordeum (no separation by soil types, only by compartment).
##### So top20 orders in rhizosphere and top20 in root for Lotus or Hordeum would be identified separately.
##### The unique orders would be combined, and one could additionally filter for RA of min. 1% in at least one of the compartments in one of the plants.

## NOTE: The Lotus nodule compartment was not taken into account, as it is mainly dominated by Rhizobiales and top20 would represent orders with very small RA.

# Step 1: Filter for Rhizosphere and Root compartments only.
Lotus_RR <- Lotus_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

Hordeum_RR <- Hordeum_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

# Step 2: Calculate mean RA per Order across all samples of each compartment (soil samples combined).
Lotus_meanRA <- Lotus_RR %>%
  group_by(Compartment, Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

Hordeum_meanRA <- Hordeum_RR %>%
  group_by(Compartment, Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

# Step 3: Identify top 20 Orders in each compartment by MeanRA.
Lotus_top20_orders_per_comp <- Lotus_meanRA %>%
  group_by(Compartment) %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20) %>%
  ungroup()

Hordeum_top20_orders_per_comp <- Hordeum_meanRA %>%
  group_by(Compartment) %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Step 4: Make a unique list of Orders across the two compartments
Lotus_unique_top20_orders <- Lotus_top20_orders_per_comp %>%
  distinct(Order) %>%
  pull(Order)

Hordeum_unique_top20_orders <- Hordeum_top20_orders_per_comp %>%
  distinct(Order) %>%
  pull(Order)

# Step 5: Check which of these Orders have a mean RA >= 0.01 in either Rhizosphere or Root.
Lotus_top_orders_above1perc <- Lotus_meanRA %>%
  filter(Order %in% Lotus_unique_top20_orders & MeanRA >= 0.01) %>%
  distinct(Order) %>%
  pull(Order)

Hordeum_top_orders_above1perc <- Lotus_meanRA %>%
  filter(Order %in% Hordeum_unique_top20_orders & MeanRA >= 0.01) %>%
  distinct(Order) %>%
  pull(Order)

# Step 6: Check how many unique orders these would be to display when combining Lotus and Hordeum top20 rhizo and/or root with >=1% RA.
combined_top_orders_1 <- unique(c(Lotus_top_orders_above1perc, Hordeum_top_orders_above1perc))
length(combined_top_orders_1) ##### 18 unique orders for display

### 3) One can go by the top20 bacterial orders across all samples, identified for Lotus and Hordeum separately.
##### One would not look for top20 separately for soil or compartment, but average it across all samples.
##### Then one could filter for orders with mean RA of min. 1%, and combine the unique orders for Lotus and Hordeum and display those.

## NOTE: Lotus nodule compartment again not taken into account, so average in rhizosphere and root calculated.
## I did test taking the Nodule compartment into account, but the result looked exactly the same, apart from loosing one group when Nodule samples were included. Therefore I kept it excluded.

# Step 1: Keep only Rhizosphere and Root compartments
Lotus_filtered <- Lotus_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

Hordeum_filtered <- Hordeum_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

# Step 2: Calculate overall mean RA per Order across these samples
Lotus_meanRA <- Lotus_filtered %>%
  group_by(Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

Hordeum_meanRA <- Hordeum_filtered %>%
  group_by(Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

# Step 3: Identify top 20 Orders by overall mean RA
Lotus_top20_orders <- Lotus_meanRA %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20)

Hordeum_top20_orders <- Hordeum_meanRA %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20)

# Step 4: Keep only those with mean RA ≥ 0.01
Lotus_top20_above1perc <- Lotus_top20_orders %>%
  filter(MeanRA >= 0.01) %>%
  pull(Order)

Hordeum_top20_above1perc <- Hordeum_top20_orders %>%
  filter(MeanRA >= 0.01) %>%
  pull(Order)

# Step 5: Combine orders Lotus and Hordeum and show how many unique orders would be displayed.
combined_top_orders_2 <- unique(c(Lotus_top20_above1perc, Hordeum_top20_above1perc))
length(combined_top_orders_2) ##### 18 unique orders for display

####### I WILL GO AHEAD WITH OPTION 3 (combined_top_orders_2).

# We will work with the selected bacterial orders, go back to the original dataframes and rename all other bacterial orders "Other".
Lotus_df.long <- Lotus_df.long %>%
  mutate(Order = if_else(Order %in% combined_top_orders_2, Order, "Other"))

Hordeum_df.long <- Hordeum_df.long %>%
  mutate(Order = if_else(Order %in% combined_top_orders_2, Order, "Other"))

# I will add a new column "Plant" to both dataframes and combine them.
Lotus_df.long <- Lotus_df.long %>%
  mutate(Plant = "Lotus")

Hordeum_df.long <- Hordeum_df.long %>%
  mutate(Plant = "Hordeum")

combined_df <- bind_rows(Lotus_df.long, Hordeum_df.long)

# Now summarise mean RA per plant - compartment - soil combination.
# df.mean_order <- combined_df %>%
#   group_by(Plant, Compartment, Soil, Order) %>%
#   summarise(RA = mean(RA, na.rm = TRUE), .groups = "drop") %>%
#   mutate(Order = factor(Order, levels = c(sort(unique(Order[Order != "Other"])), "Other")))
df.sample_order <- combined_df %>%
  group_by(sampleID, Plant, Compartment, Soil, Order) %>%
  summarise(RA=sum(RA), .groups="drop")

df.mean_order <- df.sample_order %>%
  group_by(Plant, Compartment, Soil, Order) %>%
  summarise(RA=mean(RA), .groups="drop") %>%
  mutate(Order=factor(Order, levels=c(sort(unique(Order[Order!="Other"])), "Other")),
         Plant=factor(Plant, levels=c("Lotus","Hordeum")),
         Compartment=factor(Compartment, levels=c("Rhizosphere","Root","Nodules")),
         Soil=factor(Soil, levels=c("NPK","PK","UF")))


# Set main parameters for stacked barplot.
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

# Set colours for bacterial orders.
colors <- c(
  "Acidobacteriales"   = "#570861",   # deep purple
  "Burkholderiales"    = "#645394",   # purple-blue
  "Caulobacterales"    = "#8e3563",   # magenta
  "Chitinophagales"    = "#b55385",   # rose
  "Chloroflexales"     = "#CC99BB",   # light purple-pink
  "Corynebacteriales"  = "#f6cefc",   # very light pink
  "Flavobacteriales"   = "#05294a",   # navy
  "Frankiales"         = "#114477",   # dark teal-blue
  "Gaiellales"         = "#4477AA",   # medium blue
  "Gemmatimonadales"   = "#77AADD",   # light blue
  "MB-A2-108"          = "#117777",   # teal
  "Micrococcales"      = "#44AAAA",   # turquoise,
  "Micromonosporales"  = "#99D6DD",   
  "Nitrospirales"      = "#daf0ee",   # pale aqua
  "Pedosphaerales"     = "#013220",   # very dark green
  "Propionibacteriales"= "#117744",   # forest green
  "Pseudomonadales"    = "#88CCAA",   # pastel green
  "Pseudonocardiales"  = "#95bb72",   # lime green (stays in the green cluster)
  "Rhizobiales"        = "#fdbb6b",
  "S085"               = "#774411",   # brown
  "Solibacterales"     = "#DDAA77",   # beige-brown
  "Sphingomonadales"   = "lightyellow",
  "Streptomycetales"    = "#fed5a4",   # pink-magenta (to match other actinobacteria hues)
  "Subgroup_7"         = "#AA4455",   # dark red
  "TK10"               = "#DD7788",   # reddish-pink
  "Xanthomonadales"    = "#ffc0cb",   # light pink
  "Unknown"            = "darkgrey",
  "Other"              = "lightgrey"
)

# Set order for plant and compartment in plot.
# df.mean_order <- df.mean_order %>%
#   mutate(Plant = factor(Plant, levels = c("Lotus", "Hordeum")),
#          Compartment = factor(Compartment, levels = c("Rhizosphere", "Root","Nodules")))

# Make stacked barplot.
p1 <- ggplot(df.mean_order, aes(x=Soil, y=RA, fill=Order)) +
  geom_bar(stat="identity", width=0.7) +
  scale_fill_manual(values=colors) +
  scale_y_continuous(expand=c(0,0)) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill="Bacterial order") +
  xlab("") +
  guides(fill=guide_legend(nrow=21)) +
  facet_nested(~ Plant + Compartment, scales="free_x", space="free_x") +
  theme(
    axis.text.x = element_text(size=8, color="black", angle=0, vjust=1),
    strip.text.x = element_text(size=8, face="bold"),
    legend.key.size = unit(0.25, 'cm'),
    legend.margin = margin(l = -8)
  )

p1

ggsave("Soil_WT_stackedbp_meanRA.pdf", p1, width = 10, height = 6, unit = "cm")
saveRDS(p1, file = "Soil_WT_stackedbp_meanRA.rds")
saveRDS(p1, file = "../7_final_figures/Soil_WT_stackedbp_meanRA.rds")

# Unclear what remainder is for, may not work!!!

Lotus_design <- read.table("./Lotus_data/Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt", header=TRUE, sep="\t")
Lotus_asv_table <- read.table("./Lotus_data/feature-table.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE, comment.char = "", skip = 1)
Lotus_taxonomy <- read.table("./Lotus_data/taxonomy.tsv", sep="\t", header=TRUE, fill=TRUE)

Hordeum_design <- read.table("./Hordeum_data/BarleyCSSP_Askov_reseq_metadata.txt", header=TRUE, sep="\t")
Hordeum_asv_table <- read.table("./Hordeum_data/BarleyCSSP_Askov_reseq_ASVtable_10_4.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE, comment.char = "", skip = 1)
Hordeum_taxonomy <- read.table("./Hordeum_data/Barley_Askov_Rep_10_4_taxonomy.tsv", sep="\t", header=TRUE, fill=TRUE)

# Rename columns Feature.ID to ASVid.
colnames(Lotus_taxonomy)[colnames(Lotus_taxonomy) == "Feature.ID"] <- "ASVid"
colnames(Hordeum_taxonomy)[colnames(Hordeum_taxonomy) == "Feature.ID"] <- "ASVid"

# Spike-in was used in Lotus library, this sequence/ASV will be removed from the Lotus dataset. It is not present in the Barley dataset.
Lotus_asv_table <- Lotus_asv_table[row.names(Lotus_asv_table) != "85fa8bb918a926d97659d9b64ca6fedd", ]

# Clean-up layout of taxonomy files.
Lotus_taxonomy <- Lotus_taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right") %>%
  mutate(across(Kingdom:Species, ~sub("^.{3}", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

Hordeum_taxonomy <- Hordeum_taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right") %>%
  mutate(across(Kingdom:Species, ~sub("^.{3}", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Re-order data matrices.
Lotus_design <- Lotus_design %>%
  filter(SampleID %in% colnames(Lotus_asv_table)) ## Subset design to keep samples present in ASV table.

Lotus_asv_table <- Lotus_asv_table %>% ## Subset and reorder the ASV table to match the filtered design file and only keep ASVs present in the taxonomy file.
  select(all_of(Lotus_design$SampleID)) %>%
  filter(rownames(.) %in% Lotus_taxonomy$ASVid)

Hordeum_design <- Hordeum_design %>%
  filter(Sample_ID %in% colnames(Hordeum_asv_table)) 

Hordeum_asv_table <- Hordeum_asv_table %>%
  select(all_of(Hordeum_design$SampleID)) %>%
  filter(rownames(.) %in% Hordeum_taxonomy$ASVid)

# Convert ASV reads to relative abundances and save as new dataframe.
Lotus_asv_table <- sweep(Lotus_asv_table, 2, colSums(Lotus_asv_table), "/")
Lotus_df <- as.data.frame(Lotus_asv_table) %>%
  rownames_to_column(var="ASVid") %>%
  left_join(Lotus_taxonomy %>% select(ASVid, Order), by="ASVid")

Hordeum_asv_table <- sweep(Hordeum_asv_table, 2, colSums(Hordeum_asv_table), "/")
Hordeum_df <- as.data.frame(Hordeum_asv_table) %>%
  rownames_to_column(var="ASVid") %>%
  left_join(Hordeum_taxonomy %>% select(ASVid, Order), by="ASVid")

# Reshape to long format and filter for WT samples.
Lotus_df.long <- Lotus_df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(Lotus_design %>% select(SampleID, Soil, Genotype, Compartment), by=c("sampleID"="SampleID")) %>%
  filter(Genotype=="WT")

Hordeum_df.long <- Hordeum_df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(Hordeum_design %>% select(SampleID, Plant, Soil, Genotype, Compartment), by=c("sampleID"="SampleID")) %>%
  filter(Genotype=="WT")

# Change compartment labels to be uniform in both datasets.
Lotus_df.long <- Lotus_df.long %>%
  mutate(Compartment = case_match(Compartment,
                                  "Endosphere/Rhizoplane" ~ "Root",
                                  "Rhizosphere" ~ "Rhizosphere",
                                  "Nodules" ~ "Nodules"))

Hordeum_df.long <- Hordeum_df.long %>%
  mutate(Compartment = case_match(Compartment,
                                  "rhizo" ~ "Rhizosphere",
                                  "endo"  ~ "Root"))

# Next we want to choose which bacterial orders to use for representation of the bacterial community structure in a stacked barplot.

## For this, we will first make a summary of the relative abundances by order per sample - soil - compartment group for Lotus and Hordeum.
Lotus_order_summary <- Lotus_df.long %>%
  group_by(Order, sampleID, Soil, Compartment) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

Hordeum_order_summary <- Hordeum_df.long %>%
  group_by(Order, sampleID, Soil, Compartment) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

## Next, we look at the bacterial orders with highest relative abundances in the samples in order to choose which orders to display.

## There are different ways to do this:

### 1) One can look at the top20 bacterial orders with highest mean relative abundance in the different Soil - Compartment combinations for Lotus and Hordeum separately.
##### With this option one would however end up with a large number of different bacterial orders.
##### (Even when filtering for orders with e.g. mean relative abundance of min. 1%).

### 2) One can go by the top20 bacterial orders in either the rhizosphere OR the root of Lotus or Hordeum (no separation by soil types, only by compartment).
##### So top20 orders in rhizosphere and top20 in root for Lotus or Hordeum would be identified separately.
##### The unique orders would be combined, and one could additionally filter for RA of min. 1% in at least one of the compartments in one of the plants.

## NOTE: The Lotus nodule compartment was not taken into account, as it is mainly dominated by Rhizobiales and top20 would represent orders with very small RA.

# Step 1: Filter for Rhizosphere and Root compartments only.
Lotus_RR <- Lotus_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

Hordeum_RR <- Hordeum_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

# Step 2: Calculate mean RA per Order across all samples of each compartment (soil samples combined).
Lotus_meanRA <- Lotus_RR %>%
  group_by(Compartment, Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

Hordeum_meanRA <- Hordeum_RR %>%
  group_by(Compartment, Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

# Step 3: Identify top 20 Orders in each compartment by MeanRA.
Lotus_top20_orders_per_comp <- Lotus_meanRA %>%
  group_by(Compartment) %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20) %>%
  ungroup()

Hordeum_top20_orders_per_comp <- Hordeum_meanRA %>%
  group_by(Compartment) %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Step 4: Make a unique list of Orders across the two compartments
Lotus_unique_top20_orders <- Lotus_top20_orders_per_comp %>%
  distinct(Order) %>%
  pull(Order)

Hordeum_unique_top20_orders <- Hordeum_top20_orders_per_comp %>%
  distinct(Order) %>%
  pull(Order)

# Step 5: Check which of these Orders have a mean RA >= 0.01 in either Rhizosphere or Root.
Lotus_top_orders_above1perc <- Lotus_meanRA %>%
  filter(Order %in% Lotus_unique_top20_orders & MeanRA >= 0.01) %>%
  distinct(Order) %>%
  pull(Order)

Hordeum_top_orders_above1perc <- Lotus_meanRA %>%
  filter(Order %in% Hordeum_unique_top20_orders & MeanRA >= 0.01) %>%
  distinct(Order) %>%
  pull(Order)

# Step 6: Check how many unique orders these would be to display when combining Lotus and Hordeum top20 rhizo and/or root with >=1% RA.
combined_top_orders_1 <- unique(c(Lotus_top_orders_above1perc, Hordeum_top_orders_above1perc))
length(combined_top_orders_1) ##### 18 unique orders for display

### 3) One can go by the top20 bacterial orders across all samples, identified for Lotus and Hordeum separately.
##### One would not look for top20 separately for soil or compartment, but average it across all samples.
##### Then one could filter for orders with mean RA of min. 1%, and combine the unique orders for Lotus and Hordeum and display those.

## NOTE: Lotus nodule compartment again not taken into account, so average in rhizosphere and root calculated.
## I did test taking the Nodule compartment into account, but the result looked exactly the same, apart from loosing one group when Nodule samples were included. Therefore I kept it excluded.

# Step 1: Keep only Rhizosphere and Root compartments
Lotus_filtered <- Lotus_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

Hordeum_filtered <- Hordeum_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

# Step 2: Calculate overall mean RA per Order across these samples
Lotus_meanRA <- Lotus_filtered %>%
  group_by(Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

Hordeum_meanRA <- Hordeum_filtered %>%
  group_by(Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

# Step 3: Identify top 20 Orders by overall mean RA
Lotus_top20_orders <- Lotus_meanRA %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20)

Hordeum_top20_orders <- Hordeum_meanRA %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20)

# Step 4: Keep only those with mean RA ≥ 0.01
Lotus_top20_above1perc <- Lotus_top20_orders %>%
  filter(MeanRA >= 0.01) %>%
  pull(Order)

Hordeum_top20_above1perc <- Hordeum_top20_orders %>%
  filter(MeanRA >= 0.01) %>%
  pull(Order)

# Step 5: Combine orders Lotus and Hordeum and show how many unique orders would be displayed.
combined_top_orders_2 <- unique(c(Lotus_top20_above1perc, Hordeum_top20_above1perc))
length(combined_top_orders_2) ##### 18 unique orders for display

####### I WILL GO AHEAD WITH OPTION 3 (combined_top_orders_2).

# We will work with the selected bacterial orders, go back to the original dataframes and rename all other bacterial orders "Other".
Lotus_df.long <- Lotus_df.long %>%
  mutate(Order = if_else(Order %in% combined_top_orders_2, Order, "Other"))

Hordeum_df.long <- Hordeum_df.long %>%
  mutate(Order = if_else(Order %in% combined_top_orders_2, Order, "Other"))

# I will add a new column "Plant" to both dataframes and combine them.
Lotus_df.long <- Lotus_df.long %>%
  mutate(Plant = "Lotus")

Hordeum_df.long <- Hordeum_df.long %>%
  mutate(Plant = "Hordeum")

combined_df <- bind_rows(Lotus_df.long, Hordeum_df.long)

# Now summarise mean RA per plant - compartment - soil combination.
# df.mean_order <- combined_df %>%
#   group_by(Plant, Compartment, Soil, Order) %>%
#   summarise(RA = mean(RA, na.rm = TRUE), .groups = "drop") %>%
#   mutate(Order = factor(Order, levels = c(sort(unique(Order[Order != "Other"])), "Other")))
df.sample_order <- combined_df %>%
  group_by(sampleID, Plant, Compartment, Soil, Order) %>%
  summarise(RA=sum(RA), .groups="drop")

df.mean_order <- df.sample_order %>%
  group_by(Plant, Compartment, Soil, Order) %>%
  summarise(RA=mean(RA), .groups="drop") %>%
  mutate(Order=factor(Order, levels=c(sort(unique(Order[Order!="Other"])), "Other")),
         Plant=factor(Plant, levels=c("Lotus","Hordeum")),
         Compartment=factor(Compartment, levels=c("Rhizosphere","Root","Nodules")),
         Soil=factor(Soil, levels=c("NPK","PK","UF")))


# Set main parameters for stacked barplot.
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

# Set colours for bacterial orders.
colors <- c(
  "Acidobacteriales"   = "#570861",   # deep purple
  "Burkholderiales"    = "#645394",   # purple-blue
  "Caulobacterales"    = "#8e3563",   # magenta
  "Chitinophagales"    = "#b55385",   # rose
  "Chloroflexales"     = "#CC99BB",   # light purple-pink
  "Corynebacteriales"  = "#f6cefc",   # very light pink
  "Flavobacteriales"   = "#05294a",   # navy
  "Frankiales"         = "#114477",   # dark teal-blue
  "Gaiellales"         = "#4477AA",   # medium blue
  "Gemmatimonadales"   = "#77AADD",   # light blue
  "MB-A2-108"          = "#117777",   # teal
  "Micrococcales"      = "#44AAAA",   # turquoise,
  "Micromonosporales"  = "#99D6DD",   
  "Nitrospirales"      = "#daf0ee",   # pale aqua
  "Pedosphaerales"     = "#013220",   # very dark green
  "Propionibacteriales"= "#117744",   # forest green
  "Pseudomonadales"    = "#88CCAA",   # pastel green
  "Pseudonocardiales"  = "#95bb72",   # lime green (stays in the green cluster)
  "Rhizobiales"        = "#fdbb6b",
  "S085"               = "#774411",   # brown
  "Solibacterales"     = "#DDAA77",   # beige-brown
  "Sphingomonadales"   = "lightyellow",
  "Streptomycetales"    = "#fed5a4",   # pink-magenta (to match other actinobacteria hues)
  "Subgroup_7"         = "#AA4455",   # dark red
  "TK10"               = "#DD7788",   # reddish-pink
  "Xanthomonadales"    = "#ffc0cb",   # light pink
  "Unknown"            = "darkgrey",
  "Other"              = "lightgrey"
)

# Set order for plant and compartment in plot.
# df.mean_order <- df.mean_order %>%
#   mutate(Plant = factor(Plant, levels = c("Lotus", "Hordeum")),
#          Compartment = factor(Compartment, levels = c("Rhizosphere", "Root","Nodules")))

# Make stacked barplot.
p1 <- ggplot(df.mean_order, aes(x=Soil, y=RA, fill=Order)) +
  geom_bar(stat="identity", width=0.7) +
  scale_fill_manual(values=colors) +
  scale_y_continuous(expand=c(0,0)) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill="Bacterial order") +
  xlab("") +
  guides(fill=guide_legend(nrow=21)) +
  facet_nested(~ Plant + Compartment, scales="free_x", space="free_x") +
  theme(
    axis.text.x = element_text(size=8, color="black", angle=50, hjust=1),
    strip.text.x = element_text(size=8, face="bold"),
    legend.key.size = unit(0.4, 'cm')
  )

p1

ggsave("Soil_WT_stackedbp_meanRA.pdf", p1, width = 10, height = 6, unit = "cm")
saveRDS(p1, file = "Soil_WT_stackedbp_meanRA.rds")
saveRDS(p1, file = "../7_final_figures/Soil_WT_stackedbp_meanRA.rds")

# Make heatmap that displays these data and includes info on significance (differences between soil type in each Plant-Compartment combination.)
# Prepare data for heatmap: exclude Nodules
df.heatmap <- df.sample_order %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

# Calculate mean RA per Plant x Compartment x Soil x Order
df.heatmap <- df.heatmap %>%
  group_by(Plant, Compartment, Soil, Order) %>%
  summarise(RA = mean(RA, na.rm = TRUE), .groups="drop")

# Make Order a factor: alphabetical, with "Other" last
df.heatmap <- df.heatmap %>%
  mutate(Order = factor(Order,
                        levels = c(sort(unique(Order[Order != "Other"])), "Other")),
         Plant = factor(Plant, levels=c("Lotus","Hordeum")),
         Compartment = factor(Compartment, levels=c("Rhizosphere","Root")),
         Soil = factor(Soil, levels=c("NPK","PK","UF")))

# Initialize letters list
letters_list <- list()
plants <- unique(df.heatmap$Plant)
comps  <- unique(df.heatmap$Compartment)

# Loop over Plant x Compartment x Order for significance testing
for(pl in plants){
  for(comp in comps){
    df_sub <- df.heatmap %>% filter(Plant==pl, Compartment==comp)
    for(ord in unique(df_sub$Order)){
      df_ord <- df_sub %>% filter(Order==ord)
      df_samples <- df.sample_order %>% 
        filter(Plant==pl, Compartment==comp, Order==ord)
      soil_counts <- table(df_samples$Soil)
      # Only test if at least 2 soils have >1 sample in original df.sample_order
      if(length(soil_counts) > 1 & all(soil_counts >= 2)){
        kw <- kruskal.test(RA ~ Soil, data=df_samples)
        if(!is.na(kw$p.value) & kw$p.value < 0.05){
          # Dunn post-hoc
          dunn_res <- dunnTest(RA ~ Soil, data=df_samples, method="bh")$res
          # Convert p-adj to letters
          pvals_named <- setNames(dunn_res$P.adj, gsub(" ", "", dunn_res$Comparison))
          cld <- multcompLetters(pvals_named, threshold=0.05)
          
          # Map letters to soil levels
          soil_levels <- levels(df.heatmap$Soil)
          letter_vec <- rep(NA, length(soil_levels))
          names(letter_vec) <- soil_levels
          letter_vec[names(cld$Letters)] <- cld$Letters
          
          letters_list[[paste(pl, comp, ord, sep="_")]] <- data.frame(
            Plant = pl,
            Compartment = comp,
            Order = ord,
            Soil = soil_levels,
            letter = letter_vec,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

# Combine all letters into a dataframe
df_letters <- bind_rows(letters_list)

# Join significance letters with heatmap data
df.plot <- df.heatmap %>%
  left_join(df_letters, by=c("Plant","Compartment","Order","Soil")) %>%
  mutate(
    Plant = factor(Plant, levels=c("Lotus","Hordeum")),  # only the plants you have
    Compartment = factor(Compartment, levels=c("Rhizosphere","Root")),
    Order = factor(Order, levels = c(sort(unique(Order[Order != "Other"])), "Other"))
  )

# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C"
# "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6"
# "#6A3D9A" "#FFFF99" "#B15928"

# Define breaks and colors
breaks <- c(0, 0.005, 0.052, 0.052001, 0.15999, 0.16, 0.34, 0.64)
# colors <- c("darkblue", "deepskyblue", "white","#feedb4",
            # "#F09D00", "#FFB3B2", "#ba0319", "#902121")
colors <- c("#1F78B4", "#A6CEE3", "white","#FFFF99",
            "#FF7F00", "#FB9A99", "#E31A1C", "#902121")

# Rescale breaks to 0-1 for gradientn
values <- rescale(breaks, to = c(0,1))

letter_keep <- df.plot %>% 
  group_by(Plant, Compartment, Order) %>% 
  summarise(remove = all(letter == "a")) %>%
  ungroup()

letter_keep$remove[is.na(letter_keep$remove)] <- F
Opt <- expand.grid(Plant = c("Lotus", "Hordeum"), 
                   Compartment = c("Root", "Rhizosphere"))

for(i in 1:nrow(Opt)){
  letter_keep_sub <- letter_keep[
    letter_keep$Plant == Opt$Plant[i] & 
    letter_keep$Compartment == Opt$Compartment[i],
  ]
  orders_remove_letters <- as.character(letter_keep_sub$Order[letter_keep_sub$remove])
  df.plot$letter[
    df.plot$Plant == Opt$Plant[i] & 
    df.plot$Compartment == Opt$Compartment[i] &
    df.plot$Order %in% orders_remove_letters
  ] <- NA
}

df.plot <- df.plot %>% 
  filter(!(Order %in% c("Other", "Unknown"))) %>% 
  mutate(Order = droplevels(Order))

# df.plot$letter[
#   df.plot$Plant == "Hordeum" & df.plot$Compartment == "Root" &
#   df.plot$Order == "Pseudomonadales"
# ] <- NA
# Then in ggplot:
p_heatmap <- ggplot(df.plot, aes(x=Soil, y=Order, fill=RA)) +
  geom_tile(color="grey50") +
  geom_text(aes(label=letter), na.rm=TRUE, size=3) +
  scale_fill_gradientn(
    colors = colors,
    values = values,
    limits = c(0, max(df.plot$RA, na.rm=TRUE)),
    name = "Relative abundance"
  ) +
  scale_y_discrete(
    limits = rev(levels(df.plot$Order)),
    position = "right"
  ) +
  guides(fill = guide_colorbar(
    title.position = "right",
    # title.hjust = 0.5,
    barwidth = 12.5,      # make the bar wider
    barheight = 1       # adjust thickness
  ))+
  main_theme +
  facet_nested(~ Plant + Compartment, scales="free_x", space="free_x") +
  xlab(NULL) +
  ylab("Bacterial order") +
  theme(
    axis.text.x = element_text(size=8, angle=0, vjust=1, hjust=0.5, 
                               colour = "black"),
    # axis.text.y = element_text(size=8, color="black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=8, colour = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill="grey90", color=NA),
    strip.text = element_text(size=8, face="bold"),
    legend.text = element_text(size=8, colour = "black"),
    legend.title = element_text(size=8, colour = "black"),
    legend.position = "bottom",
    plot.margin = margin(r = 10, l = 20)
  )

p_heatmap

ggsave("Soil_WT_heatmap_orders.pdf", p_heatmap, width = 12, height = 6, unit = "cm")
saveRDS(p_heatmap, file = "Soil_WT_heatmap_orders.rds")
saveRDS(p_heatmap, file = "../7_final_figures/Soil_WT_heatmap_orders.rds")
