# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages.
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggforce)
library(ggh4x)
library(ggalluvial)

# Load data.
design <- read.table("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt", header=T, sep="\t")
taxonomy <- read.table("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv", sep="\t", header=T, fill=T)
asv_table <- read.table(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  skip = 1,
  comment.char = ""
)

# Clean-up the taxonomy file layout.
taxonomy <- taxonomy %>% rename(ASVid = Feature.ID)
taxonomy <- taxonomy %>%
  separate(Taxon, into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
           sep="; ", fill="right") %>%
  mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Convert ASV reads to relative abundances and save as new dataframe.
asv_table_norm <- sweep(asv_table, 2, colSums(asv_table), "/")
df <- as.data.frame(asv_table_norm) %>%
  rownames_to_column(var="ASVid") %>%
  left_join(taxonomy %>% select(ASVid, Order), by="ASVid")

# Reshape the dataframe to a long format and filter for the soil (bulk) samples.
df.long <- df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(design %>% select(SampleID, Soil, Genotype), by=c("sampleID"="SampleID")) %>%
  filter(Genotype=="Soil")

# Summarise the relative abundance (RA) by order per sample.
df.long_order <- df.long %>%
  group_by(Order, sampleID, Soil) %>%
  summarise(RA=sum(RA, na.rm=TRUE), .groups="drop")

# Identify the top 20 bacterial orders per soil type.
top_orders <- df.long_order %>%
  group_by(Soil, Order) %>%
  summarise(MeanRA=mean(RA), .groups="drop") %>%
  group_by(Soil) %>%
  slice_max(MeanRA, n=20) %>%
  ungroup() %>%
  pull(Order) %>%
  unique()

# Group the remaining orders as "Other" and sort alphabetically.
df.long_order <- df.long_order %>%
  mutate(Order = ifelse(Order %in% top_orders, Order, "Other")) %>%
  mutate(Order = factor(Order, levels=c(sort(unique(Order[Order!="Other"])), "Other"))) %>%
  mutate(Soil = factor(Soil, levels=c("NPK","PK","UF")))

# Define the colours for the individual bacterial orders for the upcoming plot.
colors <- read.table("../../../0_files/Bacterial_order_colors.csv", header = T, sep = ",", comment.char = "")

# Set the main theme for the plot and make the stacked barplot.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  panel.border=element_rect(colour="black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(size = 6, color="black"),
  legend.text=element_text(size = 6, color="black"),
  legend.key=element_blank(),
  axis.title.y=element_text(size = 6),
  legend.position="right",
  legend.background=element_blank(),
  text=element_text(family="sans", size = 6, color="black")
)

p1 <- ggplot(df.long_order, aes(x=sampleID, y=RA, fill=Order)) +
  geom_bar(stat="identity", width=0.5) +
  facet_nested(~Soil, scales="free_x", labeller=label_parsed) +
  scale_fill_manual(values=colors$Color, breaks = colors$Order) +
  scale_y_continuous(expand=c(0,0)) +
  main_theme +
  ylab("Relative abundance") +
  labs(fill="Bacterial order") + 
  theme(axis.text.x=element_blank(),
        strip.text.x=element_text(size = 6, face="bold"),
        axis.title.x=element_blank()) +
  guides(fill=guide_legend(ncol=1))

p1

# Save the plot.
ggsave("Hordeum_bulk_order_top20_RA_stackedbp.pdf",p1, width=8, height=6, units = "cm")
saveRDS(p1, file="Hordeum_bulk_order_top20_RA_stackedbp.rds")
saveRDS(p1, file="../5_final_figure/Hordeum_bulk_order_top20_RA_stackedbp.rds")

# Now we repeat the plot but showing the mean RA of bacterial orders across the samples of the same soil type.

## Collapse non-top20 orders as "Other" using the same top_orders.
df.long2 <- df.long %>%
  mutate(Order = ifelse(Order %in% top_orders, Order, "Other"))

## Sum the RA per sample per order.
df.sample_order <- df.long2 %>%
  group_by(sampleID, Soil, Order) %>%
  summarise(RA=sum(RA), .groups="drop")

## Calculate the mean RA per soil per order.
df.mean_order <- df.sample_order %>%
  group_by(Soil, Order) %>%
  summarise(RA=mean(RA), .groups="drop") %>%
  mutate(Order=factor(Order, levels=c(sort(unique(Order[Order!="Other"])), "Other")))

## Plot.
p2 <- ggplot(df.mean_order, aes(x=Soil, y=RA, fill=Order)) +
  geom_bar(stat="identity", width=0.5) +
  scale_fill_manual(values=colors$Color, breaks = colors$Order) +
  scale_y_continuous(expand=c(0,0)) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill="Bacterial order") + 
  xlab("") +
  guides(fill=guide_legend(ncol=1))

p2

## Save the plot.
ggsave("Hordeum_bulk_order_top20_RA_mean_stackedbp.pdf", p2, width=6, height=6, units = "cm")
saveRDS(p2, file="Hordeum_bulk_order_top20_RA_mean_stackedbp.rds")
saveRDS(p2, file="../5_final_figure/Hordeum_bulk_order_top20_RA_mean_stackedbp.rds")

# Also save plot without the legend.
p2_no_legend <- p2 + theme(legend.position = "none")
p2_no_legend

ggsave("Hordeum_bulk_order_top20_RA_mean_stackedbp_no_legend.pdf", p2_no_legend, width=6, height=6, units = "cm")
saveRDS(p2_no_legend, file="Hordeum_bulk_order_top20_RA_mean_stackedbp_no_legend.rds")
saveRDS(p2_no_legend, file="../5_final_figure/Hordeum_bulk_order_top20_RA_mean_stackedbp_no_legend.rds")
