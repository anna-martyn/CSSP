# Clean up
options(warn=-1)
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages.
library(ggplot2)
library(car)
library(multcompView)
library(dplyr)

# Load nodule count data.
nod <- read.table("LotusCSSP_AskovSoils_nodule_cts.txt", header=T, sep="\t")

# Set factor order for plot.
nod$Soil_type <- factor(nod$Soil_type, levels = c("NPK","PK","UF"))

# Set the colours for the plot.
colors <- data.frame(group=c("NPK","PK","UF"),
                     colors=c("#6F944F","#B2563C","#3C7D82"))

# Significance analysis (ANOVA and Tukey HSD) and letters/asterisks details for plot.
get_tukey_letters <- function(df){
  aov_res <- aov(pink ~ Soil_type, data=df)
  tukey <- TukeyHSD(aov_res)
  
  pvals <- tukey$Soil_type[, "p adj"]
  names(pvals) <- rownames(tukey$Soil_type)
  
  letters <- multcompLetters(pvals)$Letters
  letters[levels(df$Soil_type)]
}

letters <- get_tukey_letters(nod) # Get letters

letters_df <- data.frame(
  Soil_type = levels(nod$Soil_type),
  label = letters,
  stringsAsFactors = FALSE
)

aov_res <- aov(pink ~ Soil_type, data=nod) # ANOVA p-value and asterisk
anova_p <- summary(aov_res)[[1]][["Pr(>F)"]][1]
asterisk <- case_when(
  anova_p < 0.001 ~ "***",
  anova_p < 0.01  ~ "**",
  anova_p < 0.05  ~ "*",
  TRUE ~ NA_character_
)

y_positions <- nod %>% # Y-positions for labels
  group_by(Soil_type) %>%
  summarise(y_pos = max(pink, na.rm=TRUE), .groups="drop")

letters_df <- left_join(letters_df, y_positions, by="Soil_type")

# Load plotting functions and generate plot.
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text.x = element_text(size = 8, colour = "black"),
                    axis.text.y = element_text(size = 8), colour = "black",
                    #legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

letters_df$y_pos[c(1,3)] <- c(9.5, 10.5)
p1 <- ggplot(nod, aes(x=Soil_type, y=pink, fill=Soil_type)) +
  # geom_boxplot(width=0.3, outlier.color = NA, alpha=0.7)+
  geom_boxplot(width=0.3, alpha=0.7, outlier.size = 1.5)+
  # geom_jitter(aes(), position=position_jitterdodge(jitter.width = 0.15),
  #             size=1.5, alpha=0.3)+
  scale_fill_manual(values=as.character(colors$color)) +
  geom_text(data=letters_df, aes(x=Soil_type, y=y_pos*1.1, label=label),
            inherit.aes=FALSE, size=8/.pt) +
  # annotate("text", x=2, y=max(nod$pink, na.rm=TRUE)*1.2, label=asterisk, size=10) +
  main_theme +
  ylab("Pink nodule counts/plant")+
  scale_y_continuous()+
  theme(legend.position="none", 
   # plot.title = element_text(size = 20, face="bold"),     
        # legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 8),
        # legend.text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8, colour = "black"),
        legend.key.size = unit(1,"cm"))
p1

# Save plot.
ggsave(paste("LotusWT_pink_nod.pdf", sep=""), p1, width=3, height=6, units = "cm")
saveRDS(p1, file = "LotusWT_pink_nod.rds")
saveRDS(p1, file = "../7_final_figures/LotusWT_pink_nod.rds")
