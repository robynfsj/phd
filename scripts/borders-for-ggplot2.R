# Borders for ggplot2
# ———————————————————

# This script creates the function theme_border(), which enables the plotting
# of ggplot panel.border on desired sides only. This is required in my scripts
# as I want to create bottom and left axes for each facet on my stratigraphic 
# diagrams. In addition, I wanted y-axes to have no gap between the axis and 
# the data. This is easier with panel.border than with the actual axes. An 
# alternative would be to use the lemon package, however, this requires the use
# of the function facet_rep_grid() but I need to use the tidypaleo function
# facet_abundanceh() to make my facets.

# The function has been taken from the following script:
# Cardinal, R (2013) Simple extentions to ggplot2. Available online: 
# http://egret.psychol.cam.ac.uk/statistics/R/extensions/rnc_ggplot2_border_themes_2013_01.r
# [Accessed: 13/11/20].

# I have copied it here for longevity in case the web page is ever removed.
# More details can be found at: https://egret.psychol.cam.ac.uk/statistics/R/graphs2.html

# Usage:
# Place within ggplot theme functions and select whichever sides you are after.
# panel.border = theme_border(type = c("left", "bottom", "right", "top"))




# Requirements ------------------------------------------------------------

library(grid) # for gpar




# Code duplicated from ggplot2 source for convenience ---------------------

.pt <- 1 / 0.352777778
len0_null <- function(x) {
  if (length(x) == 0)  NULL
  else                 x
}




# Generic panel border  ---------------------------------------------------
# (can set any combination of left/right/top/bottom)

theme_border <- function(
  type = c("left", "right", "bottom", "top", "none"),
  colour = "black", size = 1, linetype = 1) {
  # use with e.g.: ggplot(...) + opts( panel.border=theme_border(type=c("bottom","left")) ) + ...
  type <- match.arg(type, several.ok = TRUE)
  structure(
    list(type = type, colour = colour, size = size, linetype = linetype),
    class = c("theme_border", "element_blank", "element")
  )
}

element_grob.theme_border <- function(
  element, x = 0, y = 0, width = 1, height = 1, type = NULL, colour = NULL, 
  size = NULL, linetype = NULL, ...) {
  if (is.null(type)) type = element$type
  xlist <- c()
  ylist <- c()
  idlist <- c()
  if ("bottom" %in% type) { # bottom
    xlist <- append(xlist, c(x, x + width))
    ylist <- append(ylist, c(y, y))
    idlist <- append(idlist, c(1, 1))
  }
  if ("top" %in% type) { # top
    xlist <- append(xlist, c(x, x + width))
    ylist <- append(ylist, c(y + height, y + height))
    idlist <- append(idlist, c(2, 2))
  }
  if ("left" %in% type) { # left
    xlist <- append(xlist, c(x, x))
    ylist <- append(ylist, c(y, y + height))
    idlist <- append(idlist, c(3, 3))
  }
  if ("right" %in% type) { # right
    xlist <- append(xlist, c(x + width, x + width))
    ylist <- append(ylist, c(y, y + height))
    idlist <- append(idlist, c(4, 4))
  }
  if (length(type) == 0 || "none" %in% type) { # blank; cannot pass absence of coordinates, so pass a single point and use an invisible line
    xlist <- c(x, x)
    ylist <- c(y, y)
    idlist <- c(5, 5)
    linetype <- "blank"
  }
  gp <- gpar(lwd = len0_null(size * .pt), 
             col = colour, 
             lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), 
                     col = element$colour, 
                     lty = element$linetype)
  polylineGrob(
    x = xlist, y = ylist, id = idlist, ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}



