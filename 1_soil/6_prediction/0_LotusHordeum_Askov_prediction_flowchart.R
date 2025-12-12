# Clean up.
options(warn=-1)
rm(list=ls())

# Load required packages.
library(Gmisc, quietly = TRUE)
library(glue)
library(htmlTable)
library(grid)
library(magrittr)

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Define a custom function to create a diamond shape, and then define the individual boxes and arrows for the flowchart.
diamondGrob <- function(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar()) {
  diamond_points <- matrix(c(0.5, 0, 1, 0.5, 0.5, 1, 0, 0.5), ncol = 2, byrow = TRUE)
  polygonGrob(x = diamond_points[,1], y = diamond_points[,2], 
              default.units = "npc", gp = gp)
}

grid.newpage()

train <- boxGrob("Training\n Data",
                 x = 0.55, y = 0.9,
                 txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
                 box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

test <- boxGrob("Test Data",
                x = 0.28, y = 0.9,
                txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
                box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

mod1 <- boxGrob("Model 1: \nPredicts NPK vs. non-NPK \nTrained on NPK, PK, and UF data",
                x = 0.28, y = 0.65, width = unit(50, "mm"),
                txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
                box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

mod2 <- boxGrob("Model 2: \nPredicts PK vs. UF \nTrained on PK and UF data",
                x = 0.78, y = 0.65, width = unit(40, "mm"),
                txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
                box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

pred <- boxGrob("Prediction on Test Data: \nNPK or non-NPK?",
                x = 0.28, y = 0.45,
                txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
                box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

# Use diamondGrob as the box_fn in boxGrob.
diamond_box <- boxGrob("NPK?", box_fn = diamondGrob,
                       x = 0.28, y = 0.3,
                       txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
                       box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

pred2 <- boxGrob("Prediction on Test Data: \nPK or UF?",
                 x = 0.78, y = 0.3,
                 txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
                 box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

concl <- boxGrob("Prediction \nconcluded",
                 x = 0.78, y = 0.1,
                 txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
                 box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

concl2 <- boxGrob("Prediction \nconcluded",
                  x = 0.28, y = 0.1,
                  txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
                  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2))

No <- boxGrob("No",
              x = 0.45, y = 0.3, width = 0.05,
              txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
              box_gp = gpar(fill = "white", col = "white"))

Yes <- boxGrob("Yes",
               x = 0.28, y = 0.22, height = 0.025,
               txt_gp = gpar(col = "black", fontsize = 8, family = "helvetica"),
               box_gp = gpar(fill = "white", col = "white"))

train; test; mod1; mod2; pred; diamond_box; pred2; concl; concl2

custom_arrow <- arrow(type = "closed", length = unit(0.3, "cm"))

cc1 <- connectGrob(train, mod1, "N", "l", arrow_obj = custom_arrow)
cc2 <- connectGrob(train, mod2, "N", "l", arrow_obj = custom_arrow)
cc3 <- connectGrob(mod1, pred, "vertical", arrow_obj = custom_arrow)

cc4 <- connectGrob(test, pred, "Z", arrow_obj = custom_arrow)
x_coords <- attr(cc4, "line")$x
y_coords <- attr(cc4, "line")$y

x_coords[2] <- x_coords[2] - unit(1, "mm")
x_coords[3] <- x_coords[3] - unit(1, "mm")
x_coords[4] <- coords(pred)$left + coords(pred)$width/2
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

# Now save the flowchart.
saveRDS(grob, file = "LotusHordeum_Askov_prediction_flowchart.rds")
