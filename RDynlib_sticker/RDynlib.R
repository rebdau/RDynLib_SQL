library(ggplot2)
library(png)
library(grid)
library(hexSticker)

#' @param x x offset of the hexagon's center
#'
#' @param y y offset of the hexagon's center
#'
#' @param radius the radius (side length) of the hexagon.
#'
#' @param from_radius from where should the segment be drawn? defaults to the center
#'
#' @param to_radius to where should the segment be drawn? defaults to the radius
#'
#' @param from_angle from which angle should we draw?
#'
#' @param to_angle to which angle should we draw?
#'
#' @param fill fill color
#'
#' @param color line color
#'
#' @param size size of the line?
hex_segment2 <- function(x = 1, y = 1, radius = 1, from_radius = 0,
                         to_radius = radius, from_angle = 30, to_angle = 90,
                         fill = NA, color = NA, linewidth = 1.2) {
  from_angle <- from_angle * pi / 180
  to_angle <- to_angle * pi / 180
  coords <- data.frame(x = x + c(from_radius * cos(from_angle),
                                 to_radius * cos(from_angle),
                                 to_radius * cos(to_angle),
                                 from_radius * cos(to_angle)),
                       y = y + c(from_radius * sin(from_angle),
                                 to_radius * sin(from_angle),
                                 to_radius * sin(to_angle),
                                 from_radius * sin(to_angle))
  )
  geom_polygon(aes(x = coords$x, y = coords$y), data = coords,
               fill = fill, color = color, linewidth = linewidth)
}


library(ggplot2)
library(hexSticker)
library(ggimage)
library(png)
library(grid)

img <- readPNG("C:/Users/ament/Desktop/Github/RDynlib_sticker/test_dy.png")
img <- rasterGrob(
  img,
  width = 1.6,
  x = 0.5,
  y = 0.5,
  interpolate = TRUE
)

## Couleurs
col_dark   <- "#008000"   # vert
col_middle <- "#4b5da5"   # bleu
col_bg     <- "#FFFFFF"

font_text <- "Aller_Rg"

hex <- ggplot() +
  geom_hexagon(size = 1.2, fill = col_bg, color = NA) +
  
  ## Fond
  hex_segment2(linewidth = 0, fill = paste0(col_middle, "10"),
               from_radius = 0, to_radius = 0.9,
               from_angle = 330, to_angle = 30) +
  hex_segment2(linewidth = 0, fill = paste0(col_dark, "20"),
               from_radius = 0, to_radius = 0.9,
               from_angle = 30, to_angle = 90) +
  hex_segment2(linewidth = 0, fill = paste0(col_dark, "20"),
               from_radius = 0, to_radius = 0.9,
               from_angle = 90, to_angle = 150) +
  hex_segment2(linewidth = 0, fill = paste0(col_middle, "10"),
               from_radius = 0, to_radius = 0.9,
               from_angle = 150, to_angle = 210) +
  hex_segment2(linewidth = 0, fill = paste0(col_middle, "10"),
               from_radius = 0, to_radius = 0.9,
               from_angle = 210, to_angle = 270) +
  hex_segment2(linewidth = 0, fill = paste0(col_middle, "10"),
               from_radius = 0, to_radius = 0.9,
               from_angle = 270, to_angle = 330) +
  
  ## Bordure
  hex_segment2(linewidth = 0, fill = paste0(col_dark, "CE"),
               from_radius = 0.9, to_radius = 1,
               from_angle = 330, to_angle = 30) +
  hex_segment2(linewidth = 0, fill = paste0(col_dark, "AA"),
               from_radius = 0.9, to_radius = 1,
               from_angle = 30, to_angle = 90) +
  hex_segment2(linewidth = 0, fill = paste0(col_dark, "AA"),
               from_radius = 0.9, to_radius = 1,
               from_angle = 90, to_angle = 150) +
  hex_segment2(linewidth = 0, fill = paste0(col_dark, "CE"),
               from_radius = 0.9, to_radius = 1,
               from_angle = 150, to_angle = 210) +
  hex_segment2(linewidth = 0, fill = paste0(col_dark, "EC"),
               from_radius = 0.9, to_radius = 1,
               from_angle = 210, to_angle = 270) +
  hex_segment2(linewidth = 0, fill = paste0(col_dark, "EC"),
               from_radius = 0.9, to_radius = 1,
               from_angle = 270, to_angle = 330) +
  
  ## Logo
  geom_subview(
    subview = img,
    x = 1,
    y = 0.93,
    width = 0.95,
    height = 0.95
  ) +
  
  theme_sticker()

## Afficher le sticker
print(hex)

## Sauvegarder
save_sticker(
  hex,
  filename = "C:/Users/ament/Desktop/Github/RDynlib_sticker/RDynlib_sticker.png"
)