# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load libraries
library(ggplot2)
library(dplyr)
library(multcompView)
library(forcats)

# Load chao1 and metadata file.
alpha <- read.table("HordeumCSSP_AskovSoils_chao1.txt",
                    sep = "\t",       
                    header = TRUE,    
                    row.names = 1,    
                    check.names = FALSE)

design <- read.table("HordeumCSSP_AskovSoils_metadata.txt",
                     sep = "\t",       
                     header = TRUE,    
                     row.names = 1,    
                     check.names = FALSE)

# Combine alpha diversity (chao1) info and metadata info in new "index" file.
index <- cbind(alpha[, 1],
               design[match(row.names(alpha), row.names(design)), ])
colnames(index)[1] <- "value"

# Filter to remove soil data.
index <- subset(index, Compartment %in% c("Rhizosphere", "Root") & Genotype == "WT")

# Factor levels for Soil
l1 <- c("NPK","PK","UF")
index$Soil <- factor(index$Soil, levels=l1)

# Factor levels and rename Compartment
index$Compartment <- factor(index$Compartment,
                            levels = c("Rhizosphere","Root"))

# Define colors
colors <- data.frame(group=l1, color=c("#341C02","#A06A37","#D2B48C"))
colors <- colors[match(l1, colors$group), ]

# Initialize empty data frame for Tukey letters and ANOVA p-values
label_df <- data.frame()

# Minimum offset above points
min_offset <- 60  # adjust as needed for visibility

for(comp in levels(index$Compartment)) {
  sub_df <- index %>% filter(Compartment == comp)
  
  # Compartment y-range
  y_range_comp <- max(sub_df$value) - min(sub_df$value)
  
  # One-way ANOVA
  ano <- aov(value ~ Soil, data=sub_df)
  anova_p <- summary(ano)[[1]][["Pr(>F)"]][1]
  sig <- if(anova_p <= 0.001) "***" else if(anova_p <= 0.01) "**" else if(anova_p <= 0.05) "*" else ""
  
  # Tukey letters
  pairwise <- TukeyHSD(ano)
  tukey_letters <- multcompLetters(pairwise$Soil[,4])$Letters
  letters_df <- data.frame(Soil = names(tukey_letters),
                           Letters = tukey_letters)
  
  # Letter positions: Max + 3% of compartment range
  summary_sub <- sub_df %>%
    group_by(Soil) %>%
    summarise(MaxValue = max(value), .groups="drop") %>%
    left_join(letters_df, by="Soil") %>%
    # mutate(y_position = pmin(MaxValue + 0.03 * y_range_comp, 2000),
    #        Compartment = comp,
    #        ANOVA_sig = sig)
    mutate(y_position = pmin(MaxValue + pmax(0.03 * y_range_comp, min_offset), 2000),
           Compartment = comp,
           ANOVA_sig = sig)
  
  label_df <- bind_rows(label_df, summary_sub)
}

# Make plot.
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 20, color = "black"),
                    legend.text = element_text(size=20, color = "black"),
                    legend.key=element_blank(),
                    axis.title.y = element_text(size = 20),
                    legend.position="none",
                    strip.text = element_text(size = 20, color = "black"),
                    legend.background=element_blank(),
                    plot.title = element_text(size=20, hjust=1))

# ANOVA star positions per compartment
anova_positions <- label_df %>%
  group_by(Compartment) %>%
  summarise(y_position = pmin(max(y_position) + 0.02 * (max(MaxValue)-min(MaxValue)), 2000), .groups="drop") %>%
  left_join(label_df %>% select(Compartment, ANOVA_sig) %>% distinct(), by="Compartment") %>%
  mutate(x=2)


# Make plot with fixed y-axis
p1 <- ggplot(index, aes(x=Soil, y=value, fill=Soil)) +
  geom_boxplot(alpha=0.7, position=position_dodge(width=0.7), outlier.color=NA, width=0.3) +
  geom_jitter(position=position_jitter(width=0, height=0.17), size=3, alpha=1) +
  scale_fill_manual(values=as.character(colors$color)) +
  labs(x="", y="Chao1 index", title="Hordeum") +
  geom_text(data=label_df, aes(x=Soil, y=y_position, label=Letters), inherit.aes=FALSE, size=6) +
  geom_text(data=anova_positions, aes(x=x, y=y_position, label=ANOVA_sig), inherit.aes=FALSE, size=8) +
  facet_wrap(vars(fct_relevel(Compartment, "Rhizosphere", "Root")), scales="fixed") +
  ylim(0, 2000) +
  main_theme+
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0))

p1

# Save the plot.
ggsave("Hordeum_AskovSoils_WT_chao1_rfd.pdf", p1, width=7, height=5)
saveRDS(p1, file = "Hordeum_AskovSoils_WT_chao1_rfd.rds")

# Save ANOVA Tukey HSD output file.
write.csv(label_df, file =  "Hordeum_AskovSoils_WT_chao1_ANOVA_TukeyHSD_rfd.csv")

