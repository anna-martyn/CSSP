# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
design <- read.table("BarleyCSSP_Askov_reseq_metadata.txt", header=T, sep="\t")
asv_table <- read.table("BarleyCSSP_Askov_reseq_ASVtable_10_4.txt", sep="\t", header=T, row.names=1, check.names=F)
taxonomy <- read.table("Barley_Askov_Rep_10_4_taxonomy.txt", sep="\t", header=T, fill=T)

# Load packages
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(multcompView)

# Clean taxonomy
taxonomy <- taxonomy %>%
  separate(Taxon, into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
           sep="; ", fill="right") %>%
  mutate(across(Kingdom:Species, ~sub("^.{3}", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Order)

# Filter design for soil samples only
design <- design %>% filter(!str_detect(Description, "Lj") & Sample_ID %in% colnames(asv_table))

# Subset ASV table
asv_table <- asv_table %>%
  select(all_of(design$Sample_ID)) %>%
  filter(rownames(.) %in% taxonomy$ASVid)

# Convert ASV counts to relative abundance
asv_table_norm <- sweep(asv_table, 2, colSums(asv_table), "/")

# Merge with taxonomy and reshape to long format
df_long <- asv_table_norm %>%
  as.data.frame() %>%
  rownames_to_column("ASVid") %>%
  left_join(taxonomy, by="ASVid") %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(design %>% select(Sample_ID, Soil, Genotype), by=c("sampleID"="Sample_ID")) %>%
  filter(Genotype=="Soil")

# Summarize per-sample RA per Order
df_order_sample <- df_long %>%
  group_by(Order, sampleID, Soil) %>%
  summarise(RA=sum(RA), .groups="drop")

# Identify top 20 Orders per soil type
top_orders <- df_order_sample %>%
  group_by(Soil, Order) %>%
  summarise(MeanRA=mean(RA), .groups="drop") %>%
  group_by(Soil) %>%
  slice_max(MeanRA, n=20) %>%
  ungroup() %>%
  pull(Order) %>%
  unique()

# Keep top orders, others as "Other"
df_order_sample <- df_order_sample %>%
  mutate(Order=factor(ifelse(Order %in% top_orders, Order, "Other"),
                      levels=c(sort(top_orders), "Other")),
         Soil=factor(Soil, levels=c("NPK","PK","UF")))

# Summarize mean and standard error for plotting
df_order_summary <- df_order_sample %>%
  group_by(Order, Soil) %>%
  summarise(Mean_RA=mean(RA), SE_RA=sd(RA)/sqrt(n()), .groups="drop")

# Significance analysis.
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

write.csv(final_results, "Barley_bulk_orders_RA_ANOVATukey.csv", row.names = FALSE)

# Merge letters with summary for plotting.
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

# Create dataframe for asterisks
df_plot_asterisk <- final_results %>%
  mutate(asterisk = ifelse(P_Value < 0.05, "*", NA)) %>%
  select(Order, asterisk) %>%
  left_join(df_order_summary %>%
              group_by(Order) %>%
              summarise(max_height = max(Mean_RA + SE_RA)),
            by="Order") %>%
  mutate(y_position = max_height + 0.01) 


# Colors
colors <- c("NPK"="#341C02","PK"="#A06A37","UF"="#D2B48C")

# Plot
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major=element_line(color="gray90"),
                    panel.border=element_rect(colour="black", fill=NA, linewidth=1),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(size=20, color="black"),
                    legend.text=element_text(size=20),
                    legend.key=element_blank(),
                    axis.title.y=element_text(size=20),
                    legend.position=c(0.95,0.9),
                    legend.background=element_rect(colour="black", fill=NA),
                    text=element_text(family="sans", size=20))

p <- ggplot(df_order_summary, aes(x=Order, y=Mean_RA, fill=Soil)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.8, alpha=0.8) +
  geom_errorbar(aes(ymin=Mean_RA-SE_RA, ymax=Mean_RA+SE_RA),
                width=0.3, position=position_dodge(width=0.9), color="black") +
  geom_text(data=df_plot_letters,
            aes(x=Order, y=Mean_RA + SE_RA + 0.005, label=Letter, fill=Soil),
            position=position_dodge(width=0.9),
            inherit.aes=FALSE,
            vjust=0, size=6, family="sans") +
  geom_text(data=df_plot_asterisk,
            aes(x=Order, y=y_position, label=asterisk),
            inherit.aes=FALSE,
            vjust=0, size=10, family="sans") +   # match axis text size
  scale_fill_manual(values=colors) +
  labs(x="", y="Relative Abundance") +
  scale_y_continuous(limits=c(0,0.17), expand=c(0,0)) +
  main_theme +
  theme(axis.text.x=element_text(size=20, angle=50, hjust=1))

p

ggsave(paste("Barley_barplot_bulk_top20_RA.pdf", sep=""), p, width=12, height=6)
saveRDS(p, file = "Barley_barplot_bulk_top20_RA.rds")

###########

# Now we only plot the significant orders.
sig_orders <- final_results %>%
  filter(P_Value < 0.05) %>%
  pull(Order)

df_order_summary_sig <- df_order_summary %>%
  filter(Order %in% sig_orders)

df_plot_letters_sig <- df_plot_letters %>%
  filter(Order %in% sig_orders)

df_plot_asterisk_sig <- df_plot_asterisk %>%
  filter(Order %in% sig_orders)

# Make a graph only showing orders with significant differences across soil types.
p_sig <- ggplot(df_order_summary_sig, aes(x=Order, y=Mean_RA, fill=Soil)) +
  geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.8, alpha=0.8) +
  geom_errorbar(aes(ymin=Mean_RA-SE_RA, ymax=Mean_RA+SE_RA),
                width=0.3, position=position_dodge(width=0.9), color="black") +
  geom_text(data=df_plot_letters_sig,
            aes(x=Order, y=Mean_RA + SE_RA + 0.005, label=Letter, fill=Soil),
            position=position_dodge(width=0.9),
            inherit.aes=FALSE,
            vjust=0, size=6, family="sans") +
  geom_text(data=df_plot_asterisk_sig,
            aes(x=Order, y=y_position, label=asterisk),
            inherit.aes=FALSE,
            vjust=0, size=10, family="sans") +
  scale_fill_manual(values=colors) +
  labs(x="", y="Relative Abundance") +
  scale_y_continuous(limits=c(0,0.17), expand=c(0,0)) +
  main_theme +
  theme(axis.text.x=element_text(size=20, angle=50, hjust=1),
        legend.position = "right",
        legend.background=element_blank())

p_sig

# Save the updated plot
ggsave(paste("Barley_barplot_bulk_top20_RA_sign.pdf", sep=""), p_sig, width=8, height=6)
saveRDS(p_sig, file = "Barley_barplot_bulk_top20_RA_sign.rds")

