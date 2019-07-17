# Rev: 2017_02_19 key alignment modified, rowsep, colsep bug
# Rev: 2017_12_19 row colors fixed
generate_taxa_heatmap <- function (data.obj, taxa.levels='Genus', taxa='All', meta.info, sam.ord=NULL, data.type='P',  prev=0.1, minp=0.002,
                                   row.col.dat='Phyla', phy.no=4, sepwidth=0.01, colsep=NULL, rowsep=NULL, sepcolor='black',
                                   white='white', colFnsC=NULL, colFnsF=NULL, Rowv=T, Colv=T, dendrogram='both', margins=c(5, 15), in.grid=F,  is.labCol=T, cexCol=1, cexRow=NULL,
                                   omas=c(1, 1, 1, 8), width=12, height=6, ann='All', return.obj=FALSE, ...) {
  df <- data.obj$meta.dat
  prop <- NULL
  if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
    data.obj$abund.list[['Species']] <- data.obj$otu.tab
    rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":",
                                                         data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
  }

  for (LOI in taxa.levels) {
    cat(LOI, "\n")

    if (LOI == 'All') {
      if (taxa[1] == 'All') {
        stop("Please specify the taxa names that will be included in the heatmap!\n")
      }
      prop <- NULL
      for (LOI2 in names(data.obj$abund.list)) {
        ct <- data.obj$abund.list[[LOI2]]
        prop0 <- t(t(ct) / colSums(ct))
        prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa), , drop=FALSE])

      }
      colnames(prop) <- colnames(prop0)
      if (nrow(prop) != length(taxa)) {
        warnings('Some taxa not found in abundance lists! Please check the names!\n')
      }

    } else {
      ct <- data.obj$abund.list[[LOI]]
      prop <- t(t(ct) / colSums(ct))

      if (taxa[1] != 'All') {
        prop <- prop[taxa, , drop=FALSE]
      } else {
        prop <- prop[matrixStats::rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
      }
    }

    # Deal with zeros
    if (data.type == 'B') {
      col.scheme <- c("lightyellow", "red")
      prop[, ] <- as.numeric(prop != 0)
      breaks <- c(-0.01, 0.01, 1.1)
    }
    if (data.type == 'P'){
      col.scheme = c(white, RColorBrewer::brewer.pal(11, "Spectral"))
      ind.temp <- prop != 0
      minp <- min(prop[prop!=0])/1.1
      prop[prop==0] <- minp
      prop <- log10(prop)
      breaks <- c(log10(minp)-0.01, seq(log10(minp)+0.01, 0, len=12))
    }
    if (data.type == 'R'){
      col.scheme <- c(white, colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(ncol(prop)-1))

      prop <- t(apply(prop, 1, function(x) {
        temp <- rank(x[x!=0])
        s <- (ncol(prop) - 1) / (max(temp) - min(temp))
        temp <- 1 + (temp - min(temp)) * s
        x[x!=0] <- temp
        x
      }))
      breaks <- seq(0, ncol(prop), len=ncol(prop)+1)
    }
  }

  phy <- sapply(strsplit(rownames(prop), ";"), function(x) x[1])

  obj <- iheatmapr::main_heatmap(prop, colors=col.scheme) %>%
    iheatmapr::add_row_annotation(phy) %>%
    iheatmapr::add_col_annotation(as.data.frame(data.obj$meta.dat[,meta.info, drop=FALSE])) %>%
    iheatmapr::add_row_clustering() %>%
    iheatmapr::add_col_clustering()

  return(obj)

}

heatmap.3 <- function(x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist, hclustfun = hclust, dendrogram = c("both", "row", "column", "none"),
                      symm = FALSE, scale = c("none", "row", "column"), na.rm = TRUE,
                      revC = identical(Colv, "Rowv"), add.expr, breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none", col = "heat.colors", colsep,
                      rowsep, sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1,
                      notecol = "cyan", na.color = par("bg"), trace = c("none", "column", "row", "both"),
                      tracecol = "cyan", hline = median(breaks), vline = median(breaks),
                      linecol = tracecol, margins = c(5, 5), ColSideColors = NULL, RowSideColors = NULL,
                      side.height.fraction = 0.3, cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 +
                        1/log10(nc), labRow = NULL, labCol = NULL, key = TRUE, keysize = 1.5,
                      density.info = c("none", "histogram", "density"), denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL,
                      xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, NumColSideColors = 1,
                      NumRowSideColors = 1, KeyValueName = "Value", ...) {
  invalid <- function(x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid))) else if (is.vector(x))
        return(all(is.na(x))) else return(FALSE)
  }
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none" else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE else if (Colv == "Rowv" && !isTRUE(Rowv))
      Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% c("both",
                                                                   "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column" else dedrogram <- "none"
        warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% c("both",
                                                                   "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row" else dendrogram <- "none"
        warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  } else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  } else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  } else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  } else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    } else colInd <- rowInd
  } else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  } else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  } else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd] else rownames(x) else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd] else colnames(x) else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  } else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16 else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks) else {
                      extreme <- max(abs(x), na.rm = TRUE)
                      breaks <- seq(-extreme, extreme, length = breaks)
                    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!is.null(ColSideColors)) {
      # if (!is.matrix(ColSideColors)) stop(''ColSideColors' must be a
      # matrix')
      if (!is.character(ColSideColors) || nrow(ColSideColors) !=
          nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      # lhei <- c(lhei[1], 0.2, lhei[2])
      lhei = c(lhei[1], side.height.fraction * NumColSideColors,
               lhei[2])
    }
    if (!is.null(RowSideColors)) {
      # if (!is.matrix(RowSideColors)) stop(''RowSideColors' must be a
      # matrix')
      if (!is.character(RowSideColors) || ncol(RowSideColors) !=
          nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1),
                    lmat[, 2] + 1)
      # lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction * NumRowSideColors,
                lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!is.null(RowSideColors)) {
    if (!is.matrix(RowSideColors)) {
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      box()
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[, rowInd, drop = F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors),
             las = 2, tick = FALSE)
      }
      box()
    }
  }
  if (!is.null(ColSideColors)) {
    if (!is.matrix(ColSideColors)) {
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      box()
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop = F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1, (dim(csc)[2] - 1)),
             colnames(ColSideColors), las = 2, tick = FALSE)
      }
      box()
    }
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  } else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr),
        axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks,
        ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) {
    # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", col = na.color,
          add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, xright = csep +
                                0.5 + sepwidth[1], ytop = ncol(x) + 1, lty = 1, lwd = 1, col = sepcolor,
                              border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) -
                                0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 -
                                sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol, lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  } else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  } else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    } else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, xaxt = "n",
          yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2) else if (scale == "column")
        mtext(side = 1, "Column Z-Score", line = 2) else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    } else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    } else title("Color Key")
  } else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}
# Rev: 2017_08_23 add automatically create phylo.obj
generate_rarefy_curve <- function(data.obj, phylo.obj = NULL, grp.name,
                                  depth = NULL, npoint = 10, iter.no = 5,
                                  measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
                                  ann = "", gg.cmd = "theme(legend.justification=c(1,0), legend.position=c(1,0))",
                                  wid = 5, hei = 5) {
  cat("Create rarefaction curves!\n")
  # Rev: 2017_08_23
  if (is.null(phylo.obj)) {
    phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows = T),
                          phy_tree(data.obj$tree), tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
  }
  if (is.null(depth)) {
    depth <- min(sample_sums(phylo.obj))
    phylo.even <- rarefy_even_depth(phylo.obj, rngseed = 12345)
  } else {
    if (depth > min(sample_sums(phylo.obj))) {
      ind <- sample_sums(phylo.obj) >= depth
      cat(sum(!ind), " samples do not have sufficient number of reads!\n")
      sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ]
      data.obj <- subset_data(data.obj, ind)
    }
    phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed = 12345)
  }
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  # Rev: 2016_12_12
  if (is.character(grp)) {
    grp <- factor(grp)
  }
  if (!is.factor(grp)) {
    stop("Rarefaction curve needs a factor!\n")
  }
  res <- NULL
  incr <- depth%/%npoint
  sink("temp.txt")
  for (dep in c(10, incr * (1:npoint))) {
    x <- 0
    for (i in 1:iter.no) {
      phylo.even <- rarefy_even_depth(phylo.obj, dep, rngseed = 12345 +
                                        i)
      x <- x + estimate_richness(phylo.even, measures = measures)
    }
    res <- rbind(res, t(x[, measures, drop = F]/iter.no))
  }
  colnames(res) <- rownames(df)
  sink()
  pdf(paste0("Alpha_diversity_Rarefaction_Curve_", ann, ".pdf"), width = wid,
      height = hei)
  for (i in 1:length(measures)) {
    measure <- measures[i]
    cat("Measure: ", measure, "\n")
    res2 <- res[(0:(npoint)) * length(measures) + i, , drop = F]
    m <- t(apply(res2, 1, function(x) tapply(x, grp, mean)))
    se <- t(apply(res2, 1, function(x) tapply(x, grp, function(y) sd(y)/sqrt(length(y)))))
    uci <- m + se
    lci <- m - se
    m <- melt(m)
    uci <- melt(uci)
    lci <- melt(lci)
    res2 <- cbind(c(10, incr * (1:npoint)), m[, 2:3], uci[, 3], lci[,
                                                                    3])
    colnames(res2) <- c("Depth", "Group", "mean", "max", "min")
    res2 <- as.data.frame(res2)
    res2$Group <- factor(res2$Group, levels = levels(grp))
    # write.table(res2, paste0('Alpha_diversity_Rarefaction_', ann, '_',
    # measure, '.txt'))
    obj <- ggplot(res2, aes(x = Depth, y = mean, color = Group, group = Group)) +
      geom_errorbar(aes(ymin = min, ymax = max), alpha = 0.5, width = 0.25,
                    position = position_dodge(0.2)) + geom_line() +
                    geom_point(size = 3,shape = 21, fill = "white") + labs(y = measure)
    if (!is.null(gg.cmd)) {
      obj <- obj + eval(parse(text = gg.cmd))
    }
    print(obj)
  }
  dev.off()
}
# Rev: 2017_02_19 key alignment modified, rowsep, colsep bug Rev:
# 2017_12_19 row colors fixed
OLD.generate_taxa_heatmap <- function(data.obj, taxa.levels = "Genus", taxa = "All",
                                  meta.info, sam.ord = NULL, data.type = "P", prev = 0.1, minp = 0.002,
                                  row.col.dat = "Phyla", phy.no = 4, sepwidth = 0.01, colsep = NULL,
                                  rowsep = NULL, sepcolor = "black", white = "white", colFnsC = NULL,
                                  colFnsF = NULL, Rowv = T, Colv = T, dendrogram = "both", margins = c(5,15),
                                  in.grid = F, is.labCol = T, cexCol = 1, cexRow = NULL, omas = c(1, 1, 1, 8),
                                  width = 12, height = 6, ann = "All", return.obj = FALSE,
                                  ...) {
  colsep0 <- colsep
  rowsep0 <- rowsep
  df <- data.obj$meta.dat
  # Determine the col/rowside color
  if (is.null(colFnsC)) {
    # colFnsC <- c(colorRampPalette(c('black', 'yellow', 'red'),
    # colorRampPalette(c('black', 'green')), colorRampPalette(c('black',
    # 'blue'))))
    # https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
    colFnsC <- c(colorRampPalette(colors = RColorBrewer::brewer.pal(9, "RdBu")),
                 colorRampPalette(colors = RColorBrewer::brewer.pal(9, "PRGn")), colorRampPalette(colors = RColorBrewer::brewer.pal(9,"RdYlBu")),
                 colorRampPalette(colors = RColorBrewer::brewer.pal(9, "RdGy")),
                 colorRampPalette(colors = RColorBrewer::brewer.pal(9, "PuOr")))
  }
  if (is.null(colFnsF)) {
    # rainbow3 <- function(x) { rainbow(x + 2)[1:x] }
    jet3 <- function(x) {
      jet(x + 2)[1:x]
    }
    # colFnsF <- c(rainbow3, jet3) colFnsF <- c(colorRampPalette(colors =
    # brewer.pal(8, 'Set1')), colorRampPalette(colors = brewer.pal(7,
    # 'Set2')), colorRampPalette(colors = brewer.pal(7,'Dark2')),
    # colorRampPalette(colors = jet(8)))
    colFnsF <- function(x) {
      if (x <= 6) {
        return(RColorBrewer::brewer.pal(7, "Set2")[1:x])
      } else {
        return(colorRampPalette(colors = jet(8))(x))
      }
    }
    colFnsF <- c(colFnsF)
  }
  if ("Species" %in% taxa.levels & !("Species" %in% names(data.obj$abund.list))) {
    data.obj$abund.list[["Species"]] <- data.obj$otu.tab
    rownames(data.obj$abund.list[["Species"]]) <- paste0("OTU", rownames(data.obj$otu.tab),
                                                         ":", data.obj$otu.name[, "Phylum"], ";",
                                                         data.obj$otu.name[,"Genus"])
  }
  for (LOI in taxa.levels) {
    cat(LOI, "\n")
    if (LOI == "All") {
      if (taxa[1] == "All") {
        stop("Please specify the taxa names that will be included in the heatmap!\n")
      }
      prop <- NULL
      for (LOI2 in names(data.obj$abund.list)) {
        ct <- data.obj$abund.list[[LOI2]]
        prop0 <- t(t(ct)/colSums(ct))
        prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa),
                                  , drop = FALSE])
      }
      colnames(prop) <- colnames(prop0)
      if (nrow(prop) != length(taxa)) {
        warnings("Some taxa not found in abundance lists! Please check the names!\n")
      }
    } else {
      ct <- data.obj$abund.list[[LOI]]
      prop <- t(t(ct)/colSums(ct))
      if (taxa[1] != "All") {
        prop <- prop[taxa, , drop = FALSE]
      } else {
        prop <- prop[matrixStats::rowMaxs(prop) > minp & rowSums(prop != 0) >
                       prev * ncol(prop), , drop = FALSE]
      }
    }
    # Sort row and column
    if (is.null(sam.ord)) {
      prop <- prop[order(rownames(prop)), , drop = FALSE]
    } else {
      prop <- prop[order(rownames(prop)), sam.ord, drop = FALSE]
    }
    if (is.labCol) {
      labCol <- colnames(prop)
    } else {
      labCol <- ""
    }
    # Deal with zeros
    if (data.type == "B") {
      col.scheme <- c("lightyellow", "red")
      prop[, ] <- as.numeric(prop != 0)
      breaks <- c(-0.01, 0.01, 1.1)
    }
    if (data.type == "P") {
      # col.scheme <- c(white, colFunc(50))
      col.scheme = c(white, RColorBrewer::brewer.pal(11, "Spectral"))
      # col.scheme = c(white, colorRampPalette(c('blue', 'yellow',
      # 'red'))(11))
      ind.temp <- prop != 0
      minp <- min(prop[prop != 0])/1.1
      prop[prop == 0] <- minp
      prop <- log10(prop)
      # breaks <- c(log10(minp)-0.01, seq(log10(minp)+0.01, 0, len=51))
      breaks <- c(log10(minp) - 0.01, seq(log10(minp) + 0.01, 0,
                                          len = 12))
    }
    if (data.type == "R") {
      col.scheme <- c(white, colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(ncol(prop) -1))
      # col.scheme <- c('white', colorRampPalette(c('green', 'black',
      # 'red'))(ncol(prop)-1))
      prop <- t(apply(prop, 1, function(x) {
        temp <- rank(x[x != 0])
        s <- (ncol(prop) - 1)/(max(temp) - min(temp))
        temp <- 1 + (temp - min(temp)) * s
        x[x != 0] <- temp
        x
      }))
      breaks <- seq(0, ncol(prop), len = ncol(prop) + 1)
    }
    if (nrow(prop) > 150) {
      labRow <- ""
      cexRow <- 1
    } else {
      labRow <- rownames(prop)
      if (is.null(cexRow)) {
        cexRow <- ifelse(0.5 * 60/nrow(prop) > 1, 1, 0.5 * 60/nrow(prop))
      }
    }
    if (is.null(colsep0)) {
      if (in.grid == T) {
        colsep <- c(0:ncol(prop))
      } else {
        colsep <- c(0, ncol(prop))
      }
    } else {
      colsep <- colsep0
    }
    if (is.null(rowsep0)) {
      if (in.grid == T) {
        rowsep <- c(0:nrow(prop))
      } else {
        rowsep <- c(0, nrow(prop))
      }
    } else {
      rowsep <- rowsep0
    }
    key.list <- list()
    colsidecol <- NULL
    i <- 0
    j <- 0
    for (keyID in meta.info) {
      if (is.null(sam.ord)) {
        x <- df[, keyID]
      } else {
        x <- df[sam.ord, keyID]
      }
      if (is.factor(x) | is.character(x)) {
        x <- factor(x)
        key.list[[keyID]] <- list(breaks = levels(x),
                                  colors = (colFnsF[i + 1][[1]])(nlevels(x)),
                                  base = NA, col.na = NA, right = F, include.lowest = F)
        i <- (i + 1)%%length(colFnsF)
        colsidecol <- cbind(colsidecol, key.list[[keyID]]$colors[x])
      } else {
        key.list[[keyID]] <- makecmap(x, n = 5, colFn = colFnsC[j +
                                                                  1][[1]])
        j <- (j + 1)%%length(colFnsC)
        colsidecol <- cbind(colsidecol, cmap(x, key.list[[keyID]]))
      }
    }
    colnames(colsidecol) <- meta.info
    # Rev: 2016_09_13
    colsidecol[is.na(colsidecol)] <- "white"
    # add aunbdance key
    if (data.type == "B") {
      prop.cmap <- list(breaks = c("Absence", "Presence"),
                        colors = c("lightyellow", "red"), base = NA, col.na = NA, right = F, include.lowest = F)
      KeyName <- "Abundance"
    }
    # if (data.type == 'P'){ prop.cmap <- makecmap(prop[ind.temp], n = 5,
    # colFn = colFunc) } if (data.type == 'R'){ prop.cmap <-
    # makecmap(as.vector(prop), n = 5, colFn = colFunc) }
    if (data.type == "P") {
      KeyName <- "log(Proportion)"
    }
    if (data.type == "R") {
      KeyName <- "Rank"
    }
    # add row col key Strong assumption of the structure of 'row' names
    if (row.col.dat == "Phyla") {
      if (LOI %in% c("Class", "Order", "Family", "Genus", "Species")) {
        if (LOI == "Species") {
          phy <- sapply(strsplit(rownames(prop), ":"), function(x) x[2])
          phy <- sapply(strsplit(phy, ";"), function(x) x[1])
        } else {
          phy <- sapply(strsplit(rownames(prop), ";"), function(x) x[1])
        }
        if (sum(is.na(phy)) == 0) {
          temp <- sort(table(phy), decr = T)
          if (length(temp) > phy.no) {
            rare.phy <- names(temp)[-(1:phy.no)]
            phy[phy %in% rare.phy] <- "Other"
            phy <- factor(phy, levels = c(names(temp)[(1:phy.no)],
                                          "Other"))
          } else {
            phy <- factor(phy)
          }
          rowsidecol <- rainbow(nlevels(phy))[phy]
          rowsidecol <- rbind(rowsidecol, rowsidecol)
          rownames(rowsidecol) <- c("", "")
          phy.cmap <- list(breaks = levels(phy), colors = rainbow(nlevels(phy)),
                           base = NA, col.na = NA, right = F, include.lowest = F)
        } else {
          rowsidecol <- NULL
        }
      } else {
        rowsidecol <- NULL
      }
    } else {
      rowsidecol <- NULL
    }
    pdf(paste0("Taxa_Heatmap_", LOI, "_", ann, ".pdf"), width = width,
        height = height)
    par(oma = omas)
    # if (data.type == 'R' | data.type == 'B') { dist2 <- dist } if
    # (data.type == 'P') { # Better clustering of taxa dist2 <- function(x)
    # { as.dist((1-cor(t(x)))/2) } dist2 <- dist } # Rev: 2016_09_13 handle
    # zero sd cases if (sum(rowSds(prop) == 0) != 0 | sum(colSds(prop) ==
    # 0) != 0) { dist2 <- dist warning('Zero sd produced! Euclidean
    # distance is used instead!\n') }
    # Pearson correlation distance
    obj <- heatmap.3(prop, Rowv = Rowv, Colv = Colv, dendrogram = dendrogram,
                     scale = "none", col = col.scheme, breaks = breaks, symbreaks = F,
                     trace = "none", margins = margins, colsep = colsep, rowsep = rowsep,
                     sepcolor = sepcolor, sepwidth = c(sepwidth, sepwidth), ColSideColors = colsidecol,
                     RowSideColors = rowsidecol, cexRow = cexRow, labRow = labRow,
                     labCol = labCol, cexCol = cexCol, key = (data.type != "B"),
                     density.info = "none", symkey = F, KeyValueName = KeyName,
                     NumColSideColors = 0.5 * length(meta.info), NumRowSideColors = 0.5,
                     ...)
    par(cex = 0.75)
    par(oma = c(0, 0, 1, 0))
    if (!is.null(rowsidecol) & row.col.dat == "Phyla") {
      y.cord <- (1/(length(meta.info) + 2)) * (0:((length(meta.info) +
                                                     2) - 1))
      vkey2(phy.cmap, "Phylum", x = 0, y = -0.2, stretch = 1.2)
    }
    y.cord <- (1/(length(meta.info))) * (0:(length(meta.info) - 1))
    k <- 1
    xs <- c(0.95, 1)
    for (keyID in meta.info) {
      x <- df[, keyID]
      if (is.factor(x) | is.character(x)) {
        x <- factor(x)
        vkey2(key.list[[keyID]], keyID, x = xs[(k - 1)%%2 + 1],
              y = y.cord[k], stretch = 1.2)
      } else {
        vkey(key.list[[keyID]], keyID, x = xs[(k - 1)%%2 + 1],
             y = y.cord[k], stretch = 1.2)
      }
      k <- k + 1
    }
    if (data.type == "B") {
      vkey2(prop.cmap, KeyName, x = xs[(k - 1)%%2 + 1], y = 1, stretch = 1.2)
    }
    # if (data.type == 'P'){ vkey(prop.cmap, 'log10(Proportion)', y=1,
    # stretch=1.2) } if (data.type == 'R'){ vkey(prop.cmap, 'Rank', y=1,
    # stretch=1.2) }
    dev.off()
  }
  if (return.obj == TRUE) {
    return(obj)
  }
}
# New: 2017_03_02 General heatmap without taxa data.obj$data,
# data.obj$meta.dat Rev: 2017_11_20, Add 'LogC', taxa.as.row
generate_heatmap <- function(data.obj, meta.info, sam.ord = NULL, data.type = "P",
                             row.col.dat = NULL, sepwidth = 0.01, colsep = NULL, rowsep = NULL,
                             taxa.as.row = TRUE, colFnsC = NULL, colFnsF = NULL, Rowv = T, Colv = T,
                             dendrogram = "both", margins = c(5, 15), in.grid = F, sepcolor = "black",
                             is.labCol = T, cexCol = 1, cexRow = NULL, omas = c(1, 1, 1, 8), width = 12,
                             height = 6, ann = "All", return.obj = FALSE, pdf = TRUE, ...) {
  colsep0 <- colsep
  rowsep0 <- rowsep
  df <- data.obj$meta.dat
  # Determine the col/rowside color Determine the col/rowside color
  if (is.null(colFnsC)) {
    # colFnsC <- c(colorRampPalette(c('black', 'yellow', 'red'),
    # colorRampPalette(c('black', 'green')), colorRampPalette(c('black',
    # 'blue'))))
    # https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
    colFunsC <- c(colorRampPalette(colors = RColorBrewer::brewer.pal(9, "RdBu")),
                  colorRampPalette(colors = RColorBrewer::brewer.pal(9, "PRGn")),
                  colorRampPalette(colors = RColorBrewer::brewer.pal(9,"RdYlBu")),
                  colorRampPalette(colors = RColorBrewer::brewer.pal(9, "RdGy")),
                  colorRampPalette(colors = RColorBrewer::brewer.pal(9, "PuOr")))
  }
  if (is.null(colFnsF)) {
    colFnsF <- function(x) {
      if (x <= 6) {
        return(brewer.pal(7, "Set2")[1:x])
      } else {
        return(colorRampPalette(colors = jet(8))(x))
      }
    }
    colFnsF <- c(colFnsF)
  }
  prop <- data.obj$data
  # Sort row and column
  if (!is.null(sam.ord)) {
    prop <- prop[, sam.ord, drop = FALSE]
  }
  if (is.labCol) {
    labCol <- colnames(prop)
  } else {
    labCol <- ""
  }
  if (data.type == "LogC") {
    col.scheme <- c("white", colorpanel(9, "aliceblue", "lightblue",
                                        "blue"))
    breaks <- c(-0.1, 0.1, 2, 4, 6, 8, 10, 12, 14, 16, 18)
  }
  # Deal with zeros
  if (data.type == "B") {
    col.scheme <- c("lightyellow", "red")
    breaks <- c(-0.01, 0.01, 1.1)
  }
  if (data.type == "P") {
    # col.scheme = brewer.pal(11, 'Spectral')
    col.scheme <- bluered(11)
    breaks <- c(seq(min(prop), max(prop), len = 12))
  }
  if (data.type == "R") {
    col.scheme <- c(white, colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(ncol(prop) -
                                                                            1))
    breaks <- seq(0, ncol(prop), len = ncol(prop) + 1)
  }
  if (nrow(prop) > 1500) {
    labRow <- ""
    cexRow <- 1
  } else {
    labRow <- rownames(prop)
    if (is.null(cexRow)) {
      cexRow <- ifelse(0.5 * 60/nrow(prop) > 1, 1, 0.5 * 60/nrow(prop))
    }
  }
  if (is.null(colsep0)) {
    if (in.grid == T) {
      colsep <- c(0:ncol(prop))
    } else {
      colsep <- c(0, ncol(prop))
    }
  } else {
    colsep <- colsep0
  }
  if (is.null(rowsep0)) {
    if (in.grid == T) {
      rowsep <- c(0:nrow(prop))
    } else {
      rowsep <- c(0, nrow(prop))
    }
  } else {
    rowsep <- rowsep0
  }
  key.list <- list()
  colsidecol <- NULL
  i <- 0
  j <- 0
  for (keyID in meta.info) {
    if (is.null(sam.ord)) {
      x <- df[, keyID]
    } else {
      x <- df[sam.ord, keyID]
    }
    if (is.factor(x)) {
      key.list[[keyID]] <- list(breaks = levels(x), colors = (colFnsF[i + 1][[1]])(nlevels(x)),
                                base = NA, col.na = NA, right = F, include.lowest = F)
      i <- (i + 1)%%length(colFnsF)
      colsidecol <- cbind(colsidecol, key.list[[keyID]]$colors[x])
    } else {
      key.list[[keyID]] <- makecmap(x, n = 5, colFn = colFnsC[j +
                                                                1][[1]])
      j <- (j + 1)%%length(colFnsC)
      colsidecol <- cbind(colsidecol, cmap(x, key.list[[keyID]]))
    }
  }
  colnames(colsidecol) <- meta.info
  # Rev: 2016_09_13
  colsidecol[is.na(colsidecol)] <- "white"
  # add aunbdance key
  if (data.type == "B") {
    prop.cmap <- list(breaks = c("0", "1"),
                      colors = c("lightyellow", "red"), base = NA, col.na = NA, right = F, include.lowest = F)
    KeyName <- "Value"
  }
  if (data.type == "P") {
    KeyName <- "Value"
  }
  if (data.type == "R") {
    KeyName <- "Rank"
  }
  if (data.type == "LogC") {
    KeyName <- "Log2(count+1)"
  }
  # add row col key
  if (!is.null(row.col.dat)) {
    phy <- df[, row.col.dat]
    phy <- factor(phy)
    rowsidecol <- rainbow(nlevels(phy))[phy]
    rowsidecol <- rbind(rowsidecol, rowsidecol)
    rownames(rowsidecol) <- c("", "")
    phy.cmap <- list(breaks = levels(phy), colors = rainbow(nlevels(phy)),
                     base = NA, col.na = NA, right = F, include.lowest = F)
  } else {
    rowsidecol <- NULL
  }
  if (pdf == TRUE) {
    pdf(paste0("Heatmap_", ann, ".pdf"), width = width, height = height)
  }
  par(oma = omas)
  # Pearson correlation distance
  if (taxa.as.row) {
    prop <- prop
  } else {
    prop <- t(prop)
  }
  obj <- heatmap.3(prop, Rowv = Rowv, Colv = Colv, dendrogram = dendrogram,
                   scale = "none", col = col.scheme, breaks = breaks, symbreaks = F,
                   trace = "none", margins = margins, colsep = colsep, rowsep = rowsep,
                   sepcolor = sepcolor, sepwidth = c(sepwidth, sepwidth), ColSideColors = colsidecol,
                   RowSideColors = rowsidecol, cexRow = cexRow, labRow = labRow, labCol = labCol,
                   cexCol = cexCol, key = (data.type != "B"), density.info = "none",
                   symkey = F, KeyValueName = KeyName, NumColSideColors = 0.5 * length(meta.info),
                   NumRowSideColors = 0.5, ...)
  par(cex = 0.75)
  par(oma = c(0, 0, 1, 0))
  if (!is.null(rowsidecol)) {
    y.cord <- (1/(length(meta.info) + 2)) * (0:((length(meta.info) +
                                                   2) - 1))
    vkey2(phy.cmap, rowsidecol, x = 0, y = -0.2, stretch = 1.2)
  }
  y.cord <- (1/(length(meta.info))) * (0:(length(meta.info) - 1))
  k <- 1
  xs <- c(0.95, 1)
  for (keyID in meta.info) {
    x <- df[, keyID]
    if (is.factor(x)) {
      vkey2(key.list[[keyID]], keyID, x = xs[(k - 1)%%2 + 1], y = y.cord[k],
            stretch = 1.2)
    } else {
      vkey(key.list[[keyID]], keyID, x = xs[(k - 1)%%2 + 1], y = y.cord[k],
           stretch = 1.2)
    }
    k <- k + 1
  }
  if (data.type == "B") {
    vkey2(prop.cmap, KeyName, x = xs[(k - 1)%%2 + 1], y = 1, stretch = 1.2)
  }
  if (pdf == TRUE) {
    dev.off()
  }
  if (return.obj == TRUE) {
    return(obj)
  }
}

generate_stacked_barplots <- function(data.obj, taxa.level, grp.name){
  prop <- prop.table(data.obj$abund.list[[taxa.level]],2)
  prop.m <- melt(prop[rev(order(rowMeans(prop))),])
  prop.m$factor1 <- data.obj$meta.dat[match(prop.m$Var2, rownames(data.obj$meta.dat)), grp.name]

  aggregate <- ggplot(prop.m, aes(factor1, value, fill = Var1, key=Var1) ) +
    geom_bar(stat="identity", position="fill") +
    guides() +
    scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(prop.m$Var1))))

  by_sample <- ggplot(prop.m, aes(Var2, value, fill = Var1) ) +
    geom_bar(stat="identity") +
    guides() + facet_grid(~factor1, scales="free", space="free_x") +
    scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(prop.m$Var1))))

  list(aggregate=aggregate, by_sample=by_sample)
}



# Rev: 2016_06_13 Still need to be revised: Multiple groups, the auto
# order should look similar for differnt groups; The legends do not
# work Rev: 2017_04_28 add col.map so the color scheme can be the same
# for different data sets Rev: 2017_11_09 add label for individual bar
# plot and change ordering Rev: 2018_03_15 add par(lwd..)  Rev:
# 2018_05_31 add return col.map and aggregate the taxa not in col.map
# into 'Other' Class Rev: 2018_06_07 many colors
OLD.generate_stacked_barplot <- function(data.obj, grp.name = NULL, taxa.levels = c("Phylum", "Family", "Genus"),
                                     agg.cutoff = 0.005, border = TRUE, order.auto = FALSE,
                                     order.method = "abund", order.taxa = FALSE, cex.names = 0.5, cex.top.lab = 0.75,
                                     separate = FALSE, cex.names2 = 1, indiv = TRUE, aggre = TRUE, hei1 = 6,
                                     wid1 = 9, hei2 = 6, wid2 = 9, margin = 10, ann = "", pdf = TRUE, col.map = NULL,
                                     many.colors = FALSE, ...) {
  if (is.null(grp.name)) {
    grp <- rep(1, nrow(data.obj$meta.dat))
  } else {
    grp <- factor(data.obj$meta.dat[, grp.name])
  }
  name.list <- list()
  col.list <- list()
  lwd.o <- par("lwd")
  # rearrange the order of taxa
  if (order.taxa) {
    for (taxa.level in taxa.levels) {
      abund0 <- data.obj$abund.list[[taxa.level]]
      abund0 <- t(t(abund0)/colSums(abund0))
      data.obj$abund.list[[taxa.level]] <- data.obj$abund.list[[taxa.level]][rev(order(rowMeans(abund0))),
                                                                             ]
    }
  }
  # create color mapping
  prop.list <- list()
  for (taxa.level in taxa.levels) {
    abund0 <- data.obj$abund.list[[taxa.level]]
    abund0 <- t(t(abund0)/colSums(abund0))
    abund1 <- abund0[rowMeans(abund0) >= agg.cutoff, , drop = F]
    abund2 <- abund0[rowMeans(abund0) < agg.cutoff, , drop = F]
    prop <- rbind(abund1, Other = colSums(abund2))
    colnames(prop) <- colnames(abund0)
    # name.list[[taxa.level]] <- rownames(abund1)
    prop.list[[taxa.level]] <- prop
    # rand.col <- c(sample(rainbow(nrow(prop)*2), nrow(prop)-1), 'gray')
    # Rev: 2017_04_28
    if (is.null(col.map[[taxa.level]])) {
      if (many.colors) {
        basic <- brewer.pal(12, "Paired")
        rand.col <- c(rep_len(c(lighter(basic, -0.2), lighter(basic,
                                                              -0.1), basic, lighter(basic, 0.1), lighter(basic, 0.2)),
                              nrow(prop) - 1), "gray")
      } else {
        rand.col <- c(rep_len(brewer.pal(12, "Paired"), nrow(prop) -
                                1), "gray")
      }
      # rand.col <- c(brewer.pal(nrow(prop) - 1, 'Paired'), 'gray')
      names(rand.col) <- rownames(prop)
      col.list[[taxa.level]] <- rand.col
    } else {
      # Rev: 2018_05_30
      if (sum(!(rownames(prop) %in% names(col.map[[taxa.level]])))) {
        warning("Color mapping does not contain all the taxa! Unmapped taxa will be aggreated into \"Other\" category!\n")
        umapped.taxa <- setdiff(rownames(prop), names(col.map[[taxa.level]]))
        Other2 <- colSums(prop[umapped.taxa, , drop = FALSE])
        prop <- prop[setdiff(rownames(prop), umapped.taxa), , drop = FALSE]
        if ("Other" %in% rownames(prop)) {
          prop["Other", ] <- prop["Other", ] + Other2
        } else {
          prop <- rbind(prop, Other = Other2)
          colnames(prop) <- colnames(data.obj$abund.list[[taxa.level]])
        }
        # name.list[[taxa.level]] <- rownames(prop)
        prop.list[[taxa.level]] <- prop
      }
      col.list[[taxa.level]] <- col.map[[taxa.level]]
    }
  }
  if (indiv == TRUE) {
    if (pdf == TRUE)
      pdf(paste0("Taxa_Stacked_Barplot_Overall_Compo_", ann, ".pdf"),
          height = hei1, width = wid1)
    if (border == FALSE) {
      lty.o <- par("lty")
      par(lty = 0)
    }
    mar.o <- par(mar = par("mar") + c(0, margin, 0, 0))
    for (taxa.level in taxa.levels) {
      # abund0 <- data.obj$abund.list[[taxa.level]] abund0 <- t(t(abund0) /
      # colSums(abund0)) abund1 <- abund0[rowMeans(abund0) >= agg.cutoff, ,
      # drop=F] abund2 <- abund0[rowMeans(abund0) < agg.cutoff, , drop=F]
      # abund1 <- abund0[name.list[[taxa.level]], , drop=FALSE] abund2 <- 1 -
      # colSums(abund1) prop <- rbind(abund1, Other=abund2) colnames(prop) <-
      # colnames(abund0)
      prop <- prop.list[[taxa.level]]
      # rand.col <- c(sample(rainbow(nrow(prop)*2), nrow(prop)-1), 'gray')
      # col.list[[taxa.level]] <- rand.col
      cex.legend = ifelse(nrow(prop) > 35, 35/nrow(prop) * 0.75,
                          0.75)
      # Rev: 2016_09_13, better ordering within group based on single linkage
      # Rev: 2017_11_10, ordering rule the same across group
      if (order.auto) {
        temp <- rev(order(rowMeans(prop)))
        prop <- do.call(cbind, tapply(1:length(grp), factor(grp),
                                      function(i) {
                                        if (length(i) >= 2) {
                                          prop.sub <- prop[, i, drop = FALSE]
                                          if (order.method == "single") {
                                            ord <- hclust(vegdist(t(prop.sub)), method = "single")$order
                                          } else {
                                            prop.dummy <- prop.sub
                                            prop.dummy[prop.dummy < 0.05] <- 0
                                            if (length(temp) >= 3) {
                                              ord <- (order(-prop.dummy[temp[1], ],
                                                            -prop.dummy[temp[2],], -prop.dummy[temp[3], ]))
                                            } else {
                                              ord <- (order(-prop.dummy[temp[1], ]))
                                            }
                                          }
                                          # temp <- rev(order(rowMeans(prop.sub)))
                                          prop.sub <- prop.sub[, ord, drop = FALSE]
                                          return(prop.sub)
                                        } else {
                                          return(prop[, i, drop = FALSE])
                                        }
                                      }))
      } else {
        prop <- prop[, order(grp), drop = FALSE]
      }
      grp2 <- sort(grp)
      if (border == FALSE) {
        if (is.null(grp.name)) {
          par(lwd = 0.25)
          barplot(prop, col = col.list[[taxa.level]][rownames(prop)],
                  ylab = "Proportion", las = 2, legend.text = rownames(prop),
                  cex.names = cex.names, space = 0,
                  args.legend = list(x = "left", bty = "n", cex = cex.legend, inset = c(-0.5, 0)),
                  main = taxa.level, ...)
        } else {
          par(lwd = 0.25)
          coord <- barplot(prop, ylim = c(0, 1.05), col = col.list[[taxa.level]][rownames(prop)],
                           ylab = "Proportion", las = 2, legend.text = rownames(prop),
                           cex.names = cex.names, space = 0,
                           args.legend = list(x = "left", bty = "n", cex = cex.legend, inset = c(-0.5, 0)),
                           main = taxa.level, ...)
          coord.text <- tapply(coord, grp2, mean)
          coord.st <- tapply(coord, grp2, max)[-(length(unique(grp2)))]
          coord.end <- tapply(coord, grp2, min)[-1]
          coord.line <- (coord.st + coord.end)/2
          text.names <- names(coord.text)
          for (i in 1:length(coord.line)) {
            abline(v = coord.line[i], lty = 1, lwd = 0.5)
          }
          for (i in 1:length(coord.text)) {
            if (i%%2 == 0) {
              pos <- 3
            } else {
              pos <- 3
            }
            text(coord.text[i], 1, text.names[i], pos = pos, cex = cex.top.lab)
          }
        }
      } else {
        if (is.null(grp.name)) {
          par(lwd = 0.25)
          barplot(prop, col = col.list[[taxa.level]][rownames(prop)],
                  ylab = "Proportion", las = 2, legend.text = rownames(prop),
                  cex.names = cex.names,
                  args.legend = list(x = "left", bty = "n", cex = cex.legend, inset = c(-0.5, 0)),
                  main = taxa.level, ...)
        } else {
          par(lwd = 0.25)
          coord <- barplot(prop, ylim = c(0, 1.05), col = col.list[[taxa.level]][rownames(prop)],
                           ylab = "Proportion", las = 2, legend.text = rownames(prop),
                           cex.names = cex.names,
                           args.legend = list(x = "left", bty = "n", cex = cex.legend, inset = c(-0.5, 0)),
                           main = taxa.level, ...)
          coord.text <- tapply(coord, grp2, mean)
          coord.st <- tapply(coord, grp2, max)[-(length(unique(grp2)))]
          coord.end <- tapply(coord, grp2, min)[-1]
          coord.line <- (coord.st + coord.end)/2
          text.names <- names(coord.text)
          for (i in 1:length(coord.line)) {
            abline(v = coord.line[i], lty = 1, lwd = 0.5)
          }
          for (i in 1:length(coord.text)) {
            if (i%%2 == 0) {
              pos <- 3
            } else {
              pos <- 3
            }
            text(coord.text[i], 1, text.names[i], pos = pos, cex = cex.top.lab)
          }
        }
      }
    }
    if (border == FALSE) {
      par(lty = lty.o)
    }
    par(mar = mar.o)
    if (pdf == TRUE)
      dev.off()
  }
  # Generate averaged over stack barplot
  if (!is.null(grp.name) & aggre == TRUE) {
    if (separate != TRUE) {
      if (pdf == TRUE)
        pdf(paste0("Taxa_Stacked_Barplot_Grouped_Compo_", ann,
                   "_combine.pdf"), height = hei2, width = wid2)
      mar.o <- par(mar = par("mar") + c(0, margin/length(taxa.levels),
                                        0, 0))
      mfrow.o <- par(mfrow = c(1, length(taxa.levels)))
      cex.legend <- 10
      for (taxa.level in taxa.levels) {
        temp <- aggregate(t(prop.list[[taxa.level]]), by = list(grp),
                          mean)
        prop <- as.matrix(temp[, -1])
        rownames(prop) <- temp[, 1]
        prop <- t(prop)
        newsize <- ifelse(nrow(prop) > 35, 35/nrow(prop) * 0.75,
                          0.75)
        cex.legend <- ifelse(cex.legend < newsize, cex.legend,
                             newsize)
        # prop <- prop[, order(grp)]
        par(lwd = 0.25)
        barplot(prop, col = col.list[[taxa.level]][rownames(prop)],
                ylab = "Proportion", las = 2, cex.names = cex.names2,
                main = taxa.level, ...)
        # legend.text=rownames(prop), args.legend=list(x='left', bty='n',
        # cex=cex.legend, inset=c(-2.2, 0)))
      }
      par(mar = c(0, 0, 0, 0))
      oma.o <- par(oma = c(0, 0, 0, 0))
      for (taxa.level in taxa.levels) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        # Rev: 2016_09_10
        legend("left", legend = rev(rownames(prop.list[[taxa.level]])),
               bty = "n", fill = rev(col.list[[taxa.level]][rownames(prop.list[[taxa.level]])]),
               cex = cex.legend)
      }
      par(mar = mar.o)
      par(oma = oma.o)
      par(mfrow = mfrow.o)
      if (pdf == TRUE)
        dev.off()
    } else {
      if (pdf == TRUE)
        pdf(paste0("Taxa_Stacked_Barplot_Grouped_Compo_", ann,
                   "_Separate.pdf"), height = hei2, width = wid2)
      mar.o <- par(mar = par("mar") + c(0, margin, 0, 0))
      for (taxa.level in taxa.levels) {
        temp <- aggregate(t(prop.list[[taxa.level]]), by = list(grp),
                          mean)
        prop <- as.matrix(temp[, -1])
        rownames(prop) <- temp[, 1]
        prop <- t(prop)
        cex.legend <- ifelse(nrow(prop) > 35, 35/nrow(prop) * 0.75,
                             0.75)
        par(lwd = 0.25)
        barplot(prop, col = col.list[[taxa.level]][rownames(prop)],
                ylab = "Proportion", las = 2, cex.names = cex.names2,
                main = taxa.level, legend.text = rownames(prop),
                args.legend = list(x = "left", bty = "n", cex = cex.legend, inset = c(-0.5, 0)), ...)
      }
      par(mar = mar.o)
      if (pdf == TRUE)
        dev.off()
    }
  }
  par(lwd = lwd.o)
  return(invisible(list(col.map = col.list)))
}
# New: 2018_07_13
generate_count_heatmap <- function(data.obj, ann = "CountHeatmaps", width = 6,
                                   height = 20, taxa.level = "OTU", cutoff = 10, occurence = 2,
                                   breaks = c(-0.1, 0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 18),
                                   col = c("white", colorpanel(13, "aliceblue", "lightblue", "blue")),
                                   margins = c(10, 5), cexRow = 0.2, cexCol = 0.4, ...) {
  if (taxa.level %in% c("Species", "OTU")) {
    otu.tab <- data.obj$otu.tab
  } else {
    otu.tab <- data.obj$abund.list[[taxa.level]]
  }
  otu.tab <- otu.tab[rowSums(otu.tab) > cutoff, ]
  otu.tab <- otu.tab[rowSums(otu.tab != 0) > occurence, ]
  pdf(paste0(ann, ".", taxa.level, ".pdf"), width = width, height = height)
  heatmap.2(log2(otu.tab + 1), scale = "none", col = col, trace = "none",
            breaks = breaks, key.xlab = "log2(count + 1)", density.info = "none",
            margins = margins, cexRow = cexRow, cexCol = cexCol, ...)
  dev.off()
}
