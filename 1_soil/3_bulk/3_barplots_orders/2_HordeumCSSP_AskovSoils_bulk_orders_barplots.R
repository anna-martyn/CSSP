# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(multcompView)

# Load the input data.
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

# Clean up the taxonomy layout.
taxonomy <- taxonomy %>% rename(ASVid = Feature.ID)
taxonomy <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right") %>%
  mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Convert the ASV counts to relative abundances in the ASV table.
asv_table_norm <- sweep(asv_table, 2, colSums(asv_table), "/")

# Merge the ASV table with the taxonomy information.
# Reshape the dataframe to a long format.
# Filter to only keep the soil samples.
df_long <- asv_table_norm %>%
  as.data.frame() %>%
  rownames_to_column("ASVid") %>%
  left_join(taxonomy %>% select (ASVid, Order), by="ASVid") %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(design %>% select(SampleID, Soil, Genotype), by=c("sampleID"="SampleID")) %>%
  filter(Genotype=="Soil")

# Summarise the per-sample RA for each bacterial order.
df_order_sample <- df_long %>%
  group_by(Order, sampleID, Soil) %>%
  summarise(RA=sum(RA), .groups="drop")

# Identify the top 20 bacterial orders per soil type.
top_orders <- df_order_sample %>%
  group_by(Soil, Order) %>%
  summarise(MeanRA=mean(RA), .groups="drop") %>%
  group_by(Soil) %>%
  slice_max(MeanRA, n=20) %>%
  ungroup() %>%
  pull(Order) %>%
  unique()

# Keep the top orders, and name the others orders "Other".
df_order_sample <- df_order_sample %>%
  mutate(Order=factor(ifelse(Order %in% top_orders, Order, "Other"),
                      levels=c(sort(top_orders), "Other")),
         Soil=factor(Soil, levels=c("NPK","PK","UF")))

# Summarise the mean and standard error for plotting.
df_order_summary <- df_order_sample %>%
  group_by(Order, Soil) %>%
  summarise(Mean_RA=mean(RA), SE_RA=sd(RA)/sqrt(n()), .groups="drop")

# Perform significance analysis using ANOVA and Tukey HSD.
order_levels <- levels(df_order_sample$Order)

final_results <- data.frame(Order = character(),
                            P_Value = numeric(),
                            NPK_Letters = character(),
                            PK_Letters = character(),
                            UF_Letters = character(),
                            stringsAsFactors = FALSE)

for (order_name in order_levels) {
  if(order_name == "Other") next
  df_order <- df_order_sample %>% filter(Order == order_name)
  
  ano <- aov(RA ~ Soil, data = df_order)
  anova_results <- summary(ano)
  p_value <- anova_results[[1]]["Soil", "Pr(>F)"]
  
  pairwise <- TukeyHSD(ano)
  Tukey.levels <- pairwise[["Soil"]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters
  
  NPK_label <- ifelse("NPK" %in% names(Tukey.labels), Tukey.labels["NPK"], NA)
  PK_label <- ifelse("PK" %in% names(Tukey.labels), Tukey.labels["PK"], NA)
  UF_label <- ifelse("UF" %in% names(Tukey.labels), Tukey.labels["UF"], NA)
  
  final_results <- rbind(final_results, data.frame(
    Order = order_name,
    P_Value = p_value,
    NPK_Letters = NPK_label,
    PK_Letters = PK_label,
    UF_Letters = UF_label
  ))
}

write.csv(final_results, "Hordeum_bulk_orders_RA_ANOVATukey.csv", row.names = FALSE)

# Merge the letters with the summary file for plotting.
df_plot_letters <- df_order_summary %>%
  left_join(final_results %>%
              pivot_longer(cols=c(NPK_Letters, PK_Letters, UF_Letters),
                           names_to="Soil_letter", values_to="Letter") %>%
              mutate(Soil = case_when(
                Soil_letter=="NPK_Letters" ~ "NPK",
                Soil_letter=="PK_Letters" ~ "PK",
                Soil_letter=="UF_Letters" ~ "UF"
              )) %>%
              select(Order, Soil, Letter),
            by=c("Order","Soil"))

# Create a dataframe for the asterisks.
df_plot_asterisk <- final_results %>%
  mutate(asterisk = ifelse(P_Value < 0.05, "*", NA)) %>%
  select(Order, asterisk) %>%
  left_join(df_order_summary %>%
              group_by(Order) %>%
              summarise(max_height = max(Mean_RA + SE_RA)),
            by="Order") %>%
  mutate(y_position = max_height + 0.01) 

# Set the colours for the plot.
colors <- c("NPK"="#6F944F","PK"="#B2563C","UF"="#3C7D82")

# Set the main theme for plotting and make the plot.
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major=element_line(color="gray90"),
                    panel.border=element_rect(colour="black", fill=NA, linewidth=1),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(size = 6, color="black"),
                    legend.text=element_text(size = 6),
                    legend.key=element_blank(),
                    legend.key.size = unit(0.25, 'cm'),
                    axis.title.y=element_text(size = 6),
                    legend.margin = margin(l = -8),
                    legend.position=c(0.95,0.9),
                    text=element_text(family="sans", size = 6))

p <- ggplot(df_order_summary, aes(x=Order, y=Mean_RA, fill=Soil)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9),
           width=0.8, alpha=0.8) +
  geom_errorbar(aes(ymin=Mean_RA-SE_RA, ymax=Mean_RA+SE_RA),
                width=0.3, position=position_dodge(width=0.9), color="black") +
  geom_text(data=df_plot_letters,
            aes(x=Order, y=Mean_RA + SE_RA + 0.005, label=Letter, fill=Soil),
            position=position_dodge(width=0.9),
            inherit.aes=FALSE,
            vjust=0, size = 6/.pt, family="sans") +
  scale_fill_manual(values=colors) +
  labs(x="", y="Relative Abundance") +
  expand_limits(y = 0) + 
  scale_y_continuous(expand=c(0,0), limits = c(0,0.15), breaks = seq(0, 0.15, 0.05)) +
  main_theme +
  theme(axis.text.x=element_text(size = 6, angle=50, hjust=1))

p

# Save the plot.
ggsave(paste("Hordeum_barplot_bulk_top20_RA.pdf", sep=""), p, width=12, height=6, units = "cm")
saveRDS(p, file = "Hordeum_barplot_bulk_top20_RA.rds")
saveRDS(p, file = "../5_final_figure/Hordeum_barplot_bulk_top20_RA.rds")

###########

# Now we only plot the significant orders.
sig_orders <- final_results %>%
  filter(P_Value < 0.05 & 
    !(NPK_Letters == "a" & 
      PK_Letters == "a" & 
      UF_Letters == "a")) %>%
  pull(Order)

df_order_summary_sig <- df_order_summary %>%
  filter(Order %in% sig_orders)

df_plot_letters_sig <- df_plot_letters %>%
  filter(Order %in% sig_orders)

df_plot_asterisk_sig <- df_plot_asterisk %>%
  filter(Order %in% sig_orders)

# Make a plot only showing orders with significant differences across soil types.
p_sig <- ggplot(df_order_summary_sig, aes(x=Order, y=Mean_RA, fill=Soil)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9),
           width=0.8, alpha=0.8) +
  geom_errorbar(aes(ymin=Mean_RA-SE_RA, ymax=Mean_RA+SE_RA),
                width=0.3, position=position_dodge(width=0.9), color="black") +
  geom_text(data=df_plot_letters_sig,
            aes(x=Order, y=Mean_RA + SE_RA + 0.005, label=Letter, fill=Soil),
            position=position_dodge(width=0.9),
            inherit.aes=FALSE,
            vjust=0, size = 6/.pt, family="sans") +
  scale_fill_manual(values=colors) +
  labs(x="", y="Relative Abundance") +
  scale_y_continuous(expand=c(0,0), limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
  main_theme +
  theme(axis.text.x=element_text(size = 6, angle=50, hjust=1),
        legend.position = "right",
        legend.background=element_blank())

p_sig

# Save the updated plot
ggsave(paste("Hordeum_barplot_bulk_top20_RA_sign.pdf", sep=""), p_sig, width=8, height=6, units = "cm")
saveRDS(p_sig, file = "Hordeum_barplot_bulk_top20_RA_sign.rds")
saveRDS(p_sig, file = "../5_final_figure/Hordeum_barplot_bulk_top20_RA_sign.rds")
