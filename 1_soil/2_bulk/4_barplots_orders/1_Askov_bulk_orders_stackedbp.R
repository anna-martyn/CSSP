# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data.
design <- read.table("BarleyCSSP_Askov_reseq_metadata.txt", header=TRUE, sep="\t")
asv_table <- read.table("BarleyCSSP_Askov_reseq_ASVtable_10_4.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
taxonomy <- read.table("Barley_Askov_Rep_10_4_taxonomy.txt", sep="\t", header=TRUE, fill=TRUE)

# Load packages.
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggforce)
library(ggh4x)
library(ggalluvial)

# Clean-up taxonomy file layout.
taxonomy <- taxonomy %>%
  separate(Taxon, into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
           sep="; ", fill="right") %>%
  mutate(across(Kingdom:Species, ~sub("^.{3}", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Remove the three Lotus soil samples from the dataset and only keep samples that are in the ASV table.
design <- design %>%
  filter(!str_detect(Description, "Lj") & Sample_ID %in% colnames(asv_table))

# Subset and reorder the ASV table to match the filtered design file and only keep ASVs present in the taxonomy file.
asv_table <- asv_table %>%
  select(all_of(design$Sample_ID)) %>%
  filter(rownames(.) %in% taxonomy$ASVid)

# Convert ASV reads to relative abundances and save as new dataframe.
asv_table_norm <- sweep(asv_table, 2, colSums(asv_table), "/")
df <- as.data.frame(asv_table_norm) %>%
  rownames_to_column(var="ASVid") %>%
  left_join(taxonomy %>% select(ASVid, Order), by="ASVid")

# Reshape to long format.
df.long <- df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(design %>% select(Sample_ID, Soil, Genotype), by=c("sampleID"="Sample_ID")) %>%
  filter(Genotype=="Soil")

# Summarize relative abundance by Order per sample.
df.long_order <- df.long %>%
  group_by(Order, sampleID, Soil) %>%
  summarise(RA=sum(RA, na.rm=TRUE), .groups="drop")

# Identify top 20 orders by mean RA across all soil samples.
top20_orders <- df.long_order %>%
  group_by(Order) %>%
  summarise(MeanRA=mean(RA, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n=20) %>%
  pull(Order)

# Group remaining orders as "Other" and order alphabetically.
df.long_order <- df.long_order %>%
  mutate(Order = ifelse(Order %in% top20_orders, Order, "Other")) %>%
  mutate(Order = factor(Order, levels=c(sort(unique(Order[Order!="Other"])), "Other"))) %>%
  mutate(Soil = factor(Soil, levels=c("NPK","PK","UF")))

# Define colours for orders.
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

# Plot stacked barplot.
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

p1 <- ggplot(df.long_order, aes(x=sampleID, y=RA, fill=Order)) +
  geom_bar(stat="identity", width=0.5) +
  facet_nested(~Soil, scales="free_x", labeller=label_parsed) +
  scale_fill_manual(values=colors) +
  scale_y_continuous(expand=c(0,0)) +
  main_theme +
  ylab("Relative abundance") +
  labs(fill="Bacterial order") + 
  theme(axis.text.x=element_blank(),
        strip.text.x=element_text(size=8, face="bold"),
        axis.title.x=element_blank()) +
  guides(fill=guide_legend(nrow=21))

p1

# Save plot.
ggsave("Barley_bulk_order_top20_RA_stackedbp.pdf",
       p1, width=8, height=6, units = "cm")
saveRDS(p1, file="Barley_bulk_order_top20_RA_stackedbp.rds")

# Now make the same plot but showing mean RA of orders across samples of same soil type.

## Collapse non-top20 ASVs as "Other" using the same top20_orders.
df.long2 <- df.long %>%
  mutate(Order = ifelse(Order %in% top20_orders, Order, "Other"))

## Sum RA per sample per order.
df.sample_order <- df.long2 %>%
  group_by(sampleID, Soil, Order) %>%
  summarise(RA=sum(RA), .groups="drop")

## Compute mean RA per soil per order.
df.mean_order <- df.sample_order %>%
  group_by(Soil, Order) %>%
  summarise(RA=mean(RA), .groups="drop") %>%
  mutate(Order=factor(Order, levels=c(sort(unique(Order[Order!="Other"])), "Other")))

## Plot.
p2 <- ggplot(df.mean_order, aes(x=Soil, y=RA, fill=Order)) +
  geom_bar(stat="identity", width=0.5) +
  scale_fill_manual(values=colors) +
  scale_y_continuous(expand=c(0,0)) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill="Bacterial order") + 
  xlab("") +
  guides(fill=guide_legend(nrow=21))

p2

## Save.
ggsave("Barley_bulk_order_top20_RA_mean_stackedbp.pdf",
       p2, width=6, height=6, units = "cm")
saveRDS(p2, file="Barley_bulk_order_top20_RA_mean_stackedbp.rds")
saveRDS(p2, file="../6_final_figure/Barley_bulk_order_top20_RA_mean_stackedbp.rds")

# Now save plot without legend also.
p2_no_legend <- p2 + theme(legend.position = "none")
p2_no_legend

ggsave("Barley_bulk_order_top20_RA_mean_stackedbp_no_legend.pdf",
       p2_no_legend, width=6, height=6, units = "cm")
saveRDS(p2_no_legend, file="Barley_bulk_order_top20_RA_mean_stackedbp_no_legend.rds")
saveRDS(p2_no_legend, 
        file="../6_final_figure/Barley_bulk_order_top20_RA_mean_stackedbp_no_legend.rds")
