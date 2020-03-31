cvs <- function(data.obj, grp.name, level="Genus", bootstrap.num = 100){
  count <- t(data.obj$abund.list[[level]])
  depth <- sapply(strsplit(colnames(count), "\\."), length)
  x <- count[, !grepl("unclassified", colnames(count)) & colSums(count != 0) >= 1]
  x[x == 0] <- 0.5
  x <- x/rowSums(x)
  z <- log(x)

  demo <- data.obj$meta.dat
  y <- demo[, grp.name][match(rownames(count), rownames(demo))]

  dyn.load(system.file("cvs", "cdmm.so", package="stromboli"))
  source(system.file("cvs", "cdmm.R", package="stromboli"))

  set.seed(23)
  n <- length(y); ntrn <- 70; nrep <- 100
  pe <- numeric(nrep); pe.lasso <- numeric(nrep)

  #for (i in 1:nrep) {
 #   itrn <- sample(n, ntrn)
  #  itst <- setdiff(1:n, itrn)
  #  ans <- cv.cdmm(y[itrn], z[itrn, ], refit=TRUE)
  #  bet <- ans$bet; int <- ans$int
  #  pe[i] <- mean((y[itst] - int - z[itst, ] %*% bet)^2)
  #  ans <- cv.cdmm(y[itrn], z[itrn, ], refit=TRUE, constr=FALSE)
  #  bet.lasso <- ans$bet; int.lasso <- ans$int
  #  pe.lasso[i] <- mean((y[itst] - int.lasso - z[itst, ] %*% bet.lasso)^2)
  #  cat("Rep.", i, "done.\n")
#  }

  set.seed(43)
  p <- ncol(z);
  nboot <- bootstrap.num
  bet.bcv <- matrix(, p, nboot)
  bet.bcv.lasso <- matrix(, p, nboot)
  for (i in 1:nboot) {
    bootid <- sample(1:length(y), replace=TRUE)
    bet.bcv[, i] <- cv.cdmm(y[bootid], z[bootid, ], refit=TRUE)$bet
    bet.bcv.lasso[, i] <- cv.cdmm(y[bootid], z[bootid, ], refit=TRUE, constr=FALSE)$bet
    cat("Boot.", i, "done.\n")
  }
  stab.prob <- stab.cdmm(y, z)$prob
  stab.prob.lasso <- stab.cdmm(y, z, constr=FALSE)$prob

  bcv.prob <- rowMeans(bet.bcv != 0)
  bcv.prob.lasso <- rowMeans(bet.bcv.lasso != 0)

  isel <- bcv.prob >= 0.7
  data.frame(genus=colnames(z)[isel], bcv.prob=bcv.prob[isel], stab.prob=stab.prob[isel])

  isel.lasso <- bcv.prob.lasso >= 0.7
  data.frame(genus=colnames(z)[isel.lasso], bcv.prob=bcv.prob.lasso[isel.lasso], stab.prob=stab.prob.lasso[isel.lasso])

  bcv.sgn <- rbind(rowMeans(bet.bcv > 0), rowMeans(bet.bcv < 0))
  taxa <- matrix(unlist(strsplit(colnames(z), "\\;")), 2)
  phyla <- taxa[1, ]; genera <- taxa[2, ]

  ##plotting
  bcv.sgn.m <- as.matrix(bcv.sgn)
  rownames(bcv.sgn.m) <- c("Positive", "Negative")
  colnames(bcv.sgn.m) <- genera
  bcv.sgn.m <- data.table::data.table(bcv.sgn.m, keep.rownames = TRUE)
  bcv.sgn.m <- data.table::melt(bcv.sgn.m, id.vars = "rn")
  bcv.sgn.m$phyla <- taxa[match(bcv.sgn.m$variable, taxa)-1]
  results <- list()
  results$taxa_plot <- ggplot2::ggplot(bcv.sgn.m) +
    ggplot2::geom_bar(ggplot2::aes(x=variable, fill=rn, y=value), stat="identity", alpha = .4) +
    ggplot2::facet_wrap(~phyla, scales = "free_x")  +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
  ans <- cdmm(y, z[, isel], 0)
  (bet <- as.numeric(ans$sol))
  int <- ans$int
  fitted <- int + drop(z[, isel] %*% bet)
  ran <- range(c(y, fitted))


  results$df <- data.frame(observed=y, fitted=fitted)

  results$fitted_plot <- ggplot2::ggplot(df, ggplot2::aes(x=observed, y=fitted)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept=0, slope=1, linetype="dashed") +
    ggplot2::labs(y=paste0("Fitted ", grp.name), x = paste0("Observed ", grp.name)) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))

  return(results)

}


randomForestTest <- function(x, y, perm.no = 999, ...) {
  iris.rf <- randomForest(x = x, y = y, importance = FALSE)
  to1 <- mean(abs((as.numeric(y) - 1) - predict(iris.rf, type = "prob")[,
                                                                        2])^2)  #mean squarred error
  to2 <- mean(y != predict(iris.rf, type = "response"))  #mean prediction error
  tp <- sapply(1:perm.no, function(i) {
    if (i%%10 == 0)
      cat(".")
    y.p <- sample(y)
    iris.rf <- randomForest(x = x, y = y.p, importance = FALSE)
    t1 <- mean(abs((as.numeric(y.p) - 1) - predict(iris.rf, type = "prob")[,
                                                                           2])^2)
    t2 <- mean(y.p != predict(iris.rf, type = "response"))
    c(t1, t2)
  })
  pv1 <- (sum(tp[1, ] <= to1) + 1)/(perm.no + 1)
  pv2 <- (sum(tp[2, ] <= to2) + 1)/(perm.no + 1)
  c(pv1 = pv1, pv2 = pv2)
}
perform_rf_test <- function(data.obj, grp.name, taxa.level = "Genus", perm.no = 999,
                            prev = 0.1, minp = 0, ann = "", ...) {
  if (taxa.level == "Species") {
    if (taxa.level %in% names(data.obj$abund.list)) {
      ct <- data.obj$abund.list[[taxa.level]]
    } else {
      # Accomodate different version
      ct <- data.obj$otu.tab
      rownames(ct) <- paste0("OTU", rownames(ct), ":", data.obj$otu.name[,
                                                                         "Phylum"], ";", data.obj$otu.name[, "Genus"])
      data.obj$abund.list[["Species"]] <- ct
    }
  } else {
    ct <- data.obj$abund.list[[taxa.level]]
  }
  prop <- t(t(ct)/colSums(ct))
  prop <- prop[rowMaxs(prop) > minp & rowSums(prop != 0) > prev * ncol(prop),
               , drop = FALSE]
  prop <- t(prop)
  grp <- data.obj$meta.dat[, grp.name]
  if (!is.factor(grp))
    stop("Current random Forest test is designed to deal with binary factor data!\n")
  cat("RF test P values  ...\n")
  obj <- randomForestTest(prop, grp, perm.no, ...)
  sink(paste0("OverallAssoc_RF_test_", taxa.level, "_", ann, ".txt"))
  cat("Prediction Probability:", obj[1], "\n")
  cat("Binary Response:", obj[2], "\n")
  sink()
}
