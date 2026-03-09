# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages.
pkg <- c("dplyr", "tidyr", "tibble", "ggplot2", "FSA", "multcompView", "ggtext")
for(pk in pkg){
  library(pk, character.only = T)
}

# Load data.
design <- read.table("../1_data/LotusSC_metadata.txt", header=TRUE, sep = "\t")
asv_table <- read.table("../1_data/LotusSC_ASVtable_nospike.tsv", sep = "\t",
                        header = TRUE, row.names = 1, check.names = FALSE,
                        comment.char = "")
taxonomy <- read.table("../1_data/LjSC_taxonomy.txt", sep="\t", header=TRUE,
                       fill = TRUE)

# Modify asv table file to only keep matched ASVs.
asv_table_matched <- asv_table[grepl("Lj", rownames(asv_table)), ]

# Check whether all matched ASVs are present in taxonomy file.
missing_asvs <- setdiff(rownames(asv_table_matched), taxonomy$ASVid)
if(length(missing_asvs) == 0){
  message("All matched ASVs are present in taxonomy.")
} else {
  warning(length(missing_asvs), " matched ASVs are missing in taxonomy:")
  print(missing_asvs)
}

# Convert ASV reads to relative abundances and add taxonomic info at order level.
asv_table_norm <- sweep(asv_table_matched, 2, colSums(asv_table_matched), "/")
df <- as.data.frame(asv_table_norm) %>%
  rownames_to_column(var="ASVid") %>%
  left_join(taxonomy %>% select(ASVid, order), by="ASVid")

# Reshape to long format and add design info.
df_long <- df %>%
  pivot_longer(cols=-c(ASVid, order), names_to="SampleID", values_to="RA") %>%
  left_join(design%>% select(SampleID, Compartment, Genotype), by="SampleID")


# Select two symbionts.
df.long_symbionts <- df_long %>%
  filter(ASVid %in% c("G37_LjNodule209", "LjRoot228"))

# Set genotype order for later plots and make mutant names italic.
colors_geno <- c(
  "WT"     = "#A9C289",
  "symrk"  = "#FEDA8B",
  "ccamk"  = "#FDB366",
  "nsp1"   = "#C0E4EF",
  "nsp2"   = "#6EA6CD"
)

genotype_labels <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

# Set main theme for plot.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border=element_rect(colour="black", fill=NA, linewidth=1),
  axis.line=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(size=8, color="black"),
  strip.text=element_text(size=8, color="black"),
  legend.text=element_text(size=8),
  legend.key=element_blank(),
  axis.title.y=element_text(size=8),
  legend.position="right",
  legend.background=element_blank(),
  text=element_text(family="sans", size=8, color="black")
)

# Perform ANOVA and Tukey HSD per ASVid-compartment combination.
df_letters <- list()
df_anova <- list()

for(comp in unique(df.long_symbionts$Compartment)) {
  for(asv in unique(df.long_symbionts$ASVid)) {
    df_sub <- df.long_symbionts %>%
      filter(Compartment==comp, ASVid==asv)
    
    # Skip if fewer than 2 genotypes with data
    if(length(unique(df_sub$Genotype)) < 2) next
    
    # ANOVA
    ano <- aov(RA ~ Genotype, data=df_sub)
    p_val <- summary(ano)[[1]][["Pr(>F)"]][1]
    
    # Tukey HSD
    tukey <- TukeyHSD(ano)
    if(any(is.na(tukey$Genotype[,4]))) {
      tukey_letters <- setNames(rep(NA, length(colors_geno)), names(colors_geno))
    } else {
      tukey_letters <- multcompView::multcompLetters(tukey$Genotype[,4])$Letters
    }
    
    # Save letters
    df_letters[[paste(asv, comp, sep="_")]] <- data.frame(
      Compartment = comp,
      ASVid = asv,
      Genotype = names(tukey_letters),
      Letter = as.vector(tukey_letters),
      stringsAsFactors = FALSE
    )
    
    # Save ANOVA p-value
    df_anova[[paste(asv, comp, sep="_")]] <- data.frame(
      Compartment = comp,
      ASVid = asv,
      P_Value = p_val,
      stringsAsFactors = FALSE
    )
  }
}

df_letters <- bind_rows(df_letters)
df_anova <- bind_rows(df_anova)

# Compute y positions for letters and asterisks.
y_positions <- df.long_symbionts %>%
  group_by(Compartment, ASVid) %>%
  summarise(y_max = max(RA, na.rm=TRUE), .groups="drop")

# Merge letters and y positions.
df_plot <- df.long_symbionts %>%
  left_join(df_letters, by=c("Compartment","ASVid","Genotype")) %>%
  left_join(y_positions, by=c("Compartment","ASVid")) %>%
  mutate(y_pos = y_max + 0.05)

#Set factor levels.
df_plot$Compartment <- factor(df_plot$Compartment, levels = c("Rhizosphere", "Root", "Nodules"))
df_plot$Genotype <- factor(df_plot$Genotype, levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2"))

# Merge ANOVA asterisks.
asterisk_df <- df_anova %>%
  left_join(y_positions, by=c("Compartment","ASVid")) %>%
  mutate(
    asterisk = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE ~ NA_character_
    ),
    y_pos_ast = y_max + 0.1
  )

# Make plot.
dodge <- position_dodge(width = 0.8)

p_box <- ggplot(df_plot, aes(x=Compartment, y=RA, fill=Genotype)) +
  geom_boxplot(position=dodge, width=0.7, outlier.shape=NA, alpha=0.7) +
  geom_jitter(position=position_jitter(width=0), size=1, alpha=0.3) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8), size=1, alpha=0.3) +
  geom_text(aes(x=Compartment, y=y_pos, label=Letter, fill=Genotype), position=dodge, size=8/.pt, inherit.aes=FALSE) +
  # geom_text(data=asterisk_df, aes(x=Compartment, y=y_pos_ast, label=asterisk), inherit.aes=FALSE, size=6) +
  facet_wrap(~ASVid) +
  scale_fill_manual(values=colors_geno, labels = genotype_labels) +
  scale_color_manual(values=colors_geno) +
  labs(x="", y="Relative Abundance", fill="Genotype") +
  main_theme +
  theme(
    # axis.text.x = element_text(angle=30, hjust=1),
    legend.text = element_markdown()
  )

p_box

# Save plot.
ggsave("LotusSynCom_symbionts_RA.pdf", plot = p_box, width = 12, height = 6)
saveRDS(p_box, file = "LotusSynCom_symbionts_RA.rds")
saveRDS(p_box, file = "../../3_final_figures/LotusSynCom_symbionts_RA.rds")
