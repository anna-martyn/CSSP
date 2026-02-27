# Load packages and set colours --------------------------
pkg <- c(
  "data.table", "ggplot2", "ggh4x", "ggpubr", "gridExtra",
  "cowplot", "ggtext", "Gmisc", "glue", "htmlTable", "grid"
)
for(pk in pkg){
  library(pk, character.only = T)
}

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define colours for the genotypes, soils, and bacterial orders of interest.
cols <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

colors <- c(NPK = "#6F944F", PK = "#B2563C", UF = "#3C7D82")

order_colors <- fread("../../0_files/Bacterial_order_colors.csv")

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Flowchart ----

# Define a custom function to create a diamond shape, and then define the individual 
# boxes and arrows for the flowchart.
diamondGrob <- function(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar()) {
  diamond_points <- matrix(
    c(0.5, 0, 1, 0.5, 0.5, 1, 0, 0.5),
    ncol = 2,
    byrow = TRUE
  )
  polygonGrob(
    x = diamond_points[,1],
    y = diamond_points[,2],
    default.units = "npc",
    gp = gp
  )
}

grid.newpage()

train <- boxGrob(
  "Training\n Data",
  x = 0.55, y = 0.9,
  txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2)
)

test <- boxGrob("Test Data",
                x = 0.28, y = 0.9,
                txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
                box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

mod1 <- boxGrob("Model 1: \nPredicts NPK vs. non-NPK \nTrained on NPK, PK, and UF data",
                x = 0.28, y = 0.65, width = unit(42, "mm"),
                txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
                box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

mod2 <- boxGrob("Model 2: \nPredicts PK vs. UF \nTrained on PK and UF data",
                x = 0.78, y = 0.65, width = unit(37, "mm"),
                txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
                box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

pred <- boxGrob("Prediction on Test Data: \nNPK or non-NPK?",
                x = 0.28, y = 0.45,
                txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
                box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

# Use diamondGrob as the box_fn in boxGrob.
diamond_box <- boxGrob("NPK?", box_fn = diamondGrob,
                       x = 0.28, y = 0.3,
                       txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
                       box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

pred2 <- boxGrob("Prediction on Test Data: \nPK or UF?",
                 x = 0.78, y = 0.3,
                 txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
                 box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

concl <- boxGrob("Prediction \nconcluded",
                 x = 0.78, y = 0.1,
                 txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
                 box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

concl2 <- boxGrob("Prediction \nconcluded",
                  x = 0.28, y = 0.1,
                  txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
                  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

No <- boxGrob("No",
              x = 0.45, y = 0.3, width = 0.05,
              txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
              box_gp = gpar(fill = "white", col = "white"))

Yes <- boxGrob("Yes",
               x = 0.28, y = 0.22, height = 0.025,
               txt_gp = gpar(col = "black", fontsize = 7, family = "helvetica"),
               box_gp = gpar(fill = "white", col = "white"))

train; test; mod1; mod2; pred; diamond_box; pred2; concl; concl2

custom_arrow <- arrow(type = "closed", length = unit(0.3, "cm"))

cc1 <- connectGrob(train, mod1, "N", "l", arrow_obj = custom_arrow)
cc2 <- connectGrob(train, mod2, "N", "l", arrow_obj = custom_arrow)
cc3 <- connectGrob(mod1, pred, "vertical", arrow_obj = custom_arrow)

cc4 <- connectGrob(test, pred, "Z", arrow_obj = custom_arrow)
x_coords <- attr(cc4, "line")$x
y_coords <- attr(cc4, "line")$y

x_coords[2] <- x_coords[2] + unit(1, "mm")
x_coords[3] <- x_coords[3] + unit(1, "mm")
x_coords[4] <- coords(pred)$left + coords(pred)$width/2 + unit(1.7, "mm")
y_coords[3] <- y_coords[3] + unit(16, "mm")
y_coords[4] <- y_coords[4] + unit(16, "mm")

new_connection <- linesGrob(x = x_coords, y = y_coords, 
                            gp = gpar(fill = "black"))

cc5 <-  connectGrob(pred, diamond_box, "vertical", arrow_obj = custom_arrow)
cc6 <-  connectGrob(mod2, pred2, "vertical", arrow_obj = custom_arrow)
cc7 <-  connectGrob(diamond_box, pred2, "horizontal", arrow_obj = custom_arrow,
                    lty_gp = gpar(col = "#e81313", fill = "#e81313"))
cc8 <-  connectGrob(diamond_box, concl2, "vertical", arrow_obj = custom_arrow,
                    lty_gp = gpar(col = "#008a0e", fill = "#008a0e"))
cc9 <-  connectGrob(pred2, concl, "vertical", arrow_obj = custom_arrow)

cc1; cc2; cc3; grid.draw(new_connection); cc5; cc6; cc7; cc8; cc9
No; Yes

grob <- grid.grab()

# # Save the flowchart.
# saveRDS(grob, file = "LotusHordeum_Askov_prediction_flowchart.rds")

Res <- fread("Prediction_results.csv")
tax_summary <- fread("Pred_taxonomic_composition.csv")
Ratio_amount <- fread("Pred_accuracy_summary.csv")

# Prediction results plot ----
Res[,":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
  Obs = factor(Obs, levels = c("NPK", "PK", "UF")),
  Pred = factor(Pred, levels = c("UF", "PK", "NPK")),
  Genotype = factor(Genotype, levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2"))
)]
Res[,Prediction := Obs == Pred]
genotype_labels_legend <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

ggplot(data = Res) +
  geom_count(aes(x = Obs, y = Pred), color = "lightgrey")+
  scale_size_continuous(range=c(1.5,15)) +
  geom_jitter(data = Res[Prediction == T],
              aes(x = Obs, y = Pred, fill = Genotype),
              position = position_jitter(width = 0.35, height = 0.35, seed = 1),
              shape = 21, stroke = 0.25) +
  geom_jitter(data = Res[Prediction == F],
              aes(x = Obs, y = Pred, fill = Genotype),
              position = position_jitter(width = 0.15, height = 0.15, seed = 1),
              shape = 21, stroke = 0.25) +
  scale_shape_manual(values = c(21, 21))+
  theme_bw() +
  labs(x = "Observed", y = "Predicted")+
  guides(size = "none", fill = guide_legend(override.aes = list(size=3))) +
  facet_grid(Host ~ Compartment) +
  scale_fill_manual(values = cols, labels = genotype_labels_legend) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -8),
        legend.box.margin = margin(0, 0, 0, -15),
        strip.background = element_rect(colour = NA),
        axis.title.y = element_text(size = 7, family = "Helvetica"),
        axis.title.x = element_text(size = 7, family = "Helvetica"),
        axis.text.y = element_text(size = 7, family = "Helvetica",
                                   colour = "black"),
        axis.text.x = element_text(size = 7, family = "Helvetica",
                                   colour = "black"),
        legend.text = element_markdown(size = 7, family = "Helvetica",
                                       colour = "black"),
        legend.title = element_text(size = 7, family = "Helvetica",
                                    colour = "black"),
        strip.text = element_text(size = 7, family = "Helvetica",
                                  face = "bold"),
        legend.key.size = unit(5, "mm"))+
  NULL -> g1; g1

# Taxonomic composition plot ----
tax_summary[!(Order %in% order_colors$Order), Order:= "Other"]
tax_summary[,":="(
  Order = droplevels(factor(Order, levels = order_colors$Order)),
  Plant = factor(Plant, levels = c("Lotus", "Hordeum")),
  Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
  Soil = factor(Soil, levels = c("NPK", "PK", "UF"))
)]
tax_summary <- tax_summary[order(Order)]

ggplot(data = tax_summary, aes(x = Soil, y = RA, fill = Order))+
  geom_bar(stat = "identity", position = "stack", linewidth = 0.1) +
  facet_wrap2(vars(Plant, Compartment), strip = strip_nested(), nrow = 1)+
  scale_fill_manual(
    values = order_colors$Color,
    breaks = order_colors$Order,
    name = "Bacterial orders"
  )+
  labs(x = NULL, y = "Cumulative Mean Relative Abundance")+
  guides(fill = guide_legend(nrow = 3, title.position = "top"))+
  scale_y_continuous(expand = c(.001, .001), limits = c(0, 0.12))+
  theme_bw()+
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(0, 0, 0, -15),
    legend.margin = margin(t = -5),
    strip.background = element_rect(colour = NA),
    legend.key.size = unit(0.25, 'cm'),
    legend.key.spacing.y = unit(0, 'cm'),
    legend.justification = c(0.75, 0),
    axis.title.y = element_text(size = 7, family = "Helvetica"),
    axis.text.y = element_text(size = 7, family = "Helvetica", colour = "black"),
    axis.text.x = element_text(size = 7, family = "Helvetica", colour = "black"),
    legend.text = element_text(size = 7, family = "Helvetica"),
    legend.title = element_text(size = 7, family = "Helvetica"),
    strip.text = element_text(size = 7, family = "Helvetica", face = "bold")
  )+
  NULL -> g2; g2

#  Combining plots ----
tg <- tableGrob(Ratio_amount, theme = ttheme_default(base_size = 7), rows = NULL)

gg1 <- plot_grid(grob, tg,
                 labels = c("A", "C"),
                 rel_heights = c(0.8, 0.2),
                 ncol = 1)
gg2 <- ggarrange(g1, g2, labels = c("B", "D"), ncol = 1, heights = c(0.5, 0.5))
gg <- ggarrange(gg1, gg2, ncol = 2)

ggsave(filename = "Figure3_Askov_prediction.pdf", plot = gg,
       width = 180, height = 200, units = "mm")
