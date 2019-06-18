vkey2 <- function(map, title = NA, side = 2, stretch = 1.4, x, y, wh) {
  opar <- par(xpd = NA)
  on.exit(par(opar))
  n <- length(map$breaks) + 1
  dy <- strheight("A")
  aspect <- diff(grconvertX(1:2, from = "inches"))/diff(grconvertY(1:2,
                                                                   from = "inches"))
  dx <- dy * aspect
  if (missing(wh)) {
    wh <- 1:(n - 1)
  }
  labs <- format(map$breaks[wh])
  maxlabwidth <- max(strwidth(labs))
  if (missing(x)) {
    x <- grconvertX(1, from = "nfc") - (2 * dx)
    if (side == 4)
      x <- x - maxlabwidth - dx
  } else {
    if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
  }
  if (missing(y))
    y <- par("usr")[3] + dy
  ybord <- y + ((0:(n - 1)) * dy * stretch)
  rect(x, ybord[-n], x + dx, ybord[-1], col = map$colors, border = NA)
  if (side == 4) {
    xtext <- x + dx
    text(x = x, y = ybord[n] + (1.5 * dy), title, adj = c(0, 0))
  }
  if (side == 2) {
    xtext <- x
    text(x = x + dx, y = ybord[n] + (1.5 * dy), title, adj = c(1, 0))
  }
  text(x = xtext, y = ybord[wh] + 0.5 * dy, labels = labs, pos = side)
}
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel ncol: Number of columns of plots nrow: Number of rows
    # needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols,
                     nrow = ceiling(numPlots/cols), byrow = TRUE)
  }
  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
###################################### New: 2018_03_09
set_colour_theme <- function(type = "stata", pal = palette(RColorBrewer::brewer.pal(5,
                                                                                    "Set1"))) {
  library(ggthemes)  # Load
  if (type == "stata") {
    scale_colour_discrete <<- function(...) scale_color_stata()
    scale_fill_discrete <<- function(...) scale_fill_stata()
  }
  if (type == "economist") {
    scale_colour_discrete <<- function(...) scale_color_economist()
    scale_fill_discrete <<- function(...) scale_fill_economist()
  }
  if (type == "wsj") {
    scale_colour_discrete <<- function(...) scale_color_wsj()
    scale_fill_discrete <- function(...) scale_fill_wsj()
  }
  if (type == "hc") {
    scale_colour_discrete <<- function(...) scale_color_hc()
    scale_fill_discrete <<- function(...) scale_fill_hc()
  }
  if (type == "default") {
    scale_fill_discrete <<- ggplot2::scale_fill_discrete
    scale_colour_discrete <<- ggplot2::scale_colour_discrete
  }
  if (type == "cb") {
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7")
    scale_fill_discrete <<- function(...) scale_fill_manual(values = cbPalette)
    scale_colour_discrete <<- function(...) scale_colour_manual(values = cbPalette)
  }
  if (type == "cbb") {
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                    "#0072B2", "#D55E00", "#CC79A7")
    scale_fill_discrete <<- function(...) scale_fill_manual(values = cbbPalette)
    scale_colour_discrete <<- function(...) scale_colour_manual(values = cbbPalette)
  }
  if (type == "customize") {
    scale_fill_discrete <<- function(...) scale_fill_manual(values = pal)
    scale_colour_discrete <<- function(...) scale_colour_manual(values = pal)
  }
}
col.func <- colorRampPalette(c("blue", "cyan", "red", "yellow"))
lighter <- function(color, factor = 0.2) {
  ## converts color to hsv, multiplies v by factor, returns colors as
  ## hexcode
  x = rgb2hsv(col2rgb(color))
  v = pmax(pmin(x[3, ] + factor, 1), 0)
  hsv(h = x[1, ], s = x[2, ], v = v)
}
