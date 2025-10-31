# Load data and packages -------------------------------------------------------
pkg <- c("data.table", "ggplot2", "scales", "grid", "RColorBrewer",
         "car", "multcompView")

for(pk in pkg){
  library(pk, character.only = T)
}

# directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load data
weight <- fread("LotusHordeum_AskovSoils_shootfw.txt")

# Set genotype colours
colors <- data.frame(
  group = c("WT", "symrk","ccamk","nsp1", "nsp2"), 
  color = c("#A9C289", "#FEDA8B", "#FDB366", "#C0E4EF", "#6EA6CD")
)

# Adjusting data ---------------------------------------------------------------
# Set correct levels
weight[,":="(
  Soil_type = factor(Soil_type, levels = c("NPK", "PK", "UF")),
  Genotype = factor(Genotype, levels = colors$group)
)]

# Split data into barley and lotus 
Barley_weight <- weight[Plant_species == "Hordeum"]
Lotus_weight <- weight[Plant_species == "Lotus"]

# Settings ---------------------------------------------------------------------
dodge <- position_dodge(width = 0.9)

main_theme <- theme(
  panel.background=element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text.x = element_text(size = 8, angle = 45, vjust = 1, 
                             hjust=1, colour = "black"),
  axis.text.y = element_text(size = 8, colour = "black"),
  #legend.position="top",
  legend.background=element_blank(),
  legend.key=element_blank(),
  text=element_text(family="sans")
)



# Matrix grid
Opt <- expand.grid(Plant = unique(weight$Plant_species),
                   Soil = unique(weight$Soil_type))

Labels_list <- as.list(rep(NA, nrow(Opt)))
for(i in 1:nrow(Opt)){
  weight_sub <- weight[
    Plant_species == Opt$Plant[i] & Soil_type == Opt$Soil[i]
  ]
  ano <- aov(Fresh_weight ~ Genotype, data = weight_sub)
  pairwise <- TukeyHSD(ano)
  
  LABELS <- multcompLetters(pairwise$Genotype[,"p adj"])
  dt <- data.table(Genotype = names(LABELS$Letters),
                   Letter = LABELS$Letters,
                   Plant_species = Opt$Plant[i],
                   Soil_type = Opt$Soil[i] )
  Labels_list[[i]] <- dt
}

Labels_dt <- rbindlist(Labels_list)

W <- weight[,.(y_pos = max(Fresh_weight)), 
            list(Plant_species, Soil_type, Genotype)]

Labels_dt <- merge(Labels_dt, W)
Labels_dt[,y_pos:=y_pos + ifelse(Plant_species == "Lotus", 0.01, 0.2)]

p1 <- ggplot(Lotus_weight, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  #geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  geom_text(data=Labels_dt[Plant_species == "Lotus"], 
            aes(x=Genotype, y = y_pos, label = Letter),
            inherit.aes=FALSE, size=8/.pt) +
  geom_boxplot(width=0.5,position = dodge,outlier.color = NA)+
  geom_jitter(position=position_jitterdodge(jitter.width = 0.3),
              size=0.5, alpha=0.3)+
  scale_fill_manual(values=as.character(colors$color)) +
  facet_wrap(~Soil_type, scales = "free_x", nrow=1)+
  main_theme +
  ylab("Shoot fresh weight/plant [g]")+
  scale_y_continuous(limits = c(0, 0.15))+
  ggtitle("Lotus")+
  theme(
    legend.position="none", 
    plot.title = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8, colour = "black"),
    strip.text = element_text(size = 8, colour = "black", face = "bold"),
    legend.text = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(size = 8, angle = 45, vjust = 1,
                               hjust=1, colour = "black"),
    legend.key.size = unit(0.5,"cm")
  )
p1

ggsave(paste("Lotus_CSSP_fw_BOX_v2.png", sep=""), 
       p1, width=6, height=5, units = "cm")

p2 <- ggplot(Barley_weight, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  #geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  geom_text(data=Labels_dt[Plant_species == "Hordeum"], 
            aes(x=Genotype, y = y_pos, label = Letter),
            inherit.aes=FALSE, size=8/.pt) +
  geom_boxplot(width=0.5,position = dodge,outlier.color = NA)+
  geom_jitter(position=position_jitterdodge(jitter.width = 0.3),
              size=2, alpha=0.3)+
  scale_fill_manual(values=as.character(colors$color)) +
  facet_wrap(~Soil_type, scales = "free_x", nrow=1)+
  main_theme +
  ylab("Shoot fresh weight/plant [g]")+
  scale_y_continuous()+
  ggtitle("Barley")+
  theme(
    legend.position="none", 
    plot.title = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8, colour = "black"),
    strip.text = element_text(size = 8, colour = "black", face = "bold"),
    legend.text = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(size = 8, angle = 45, vjust = 1,
                               hjust=1, colour = "black"),
    legend.key.size = unit(0.5,"cm")
  )
p2

ggsave(paste("Barley_CSSP_fw_BOX_v2.png", sep=""), p2,
       width=6, height=5, units = "cm")
