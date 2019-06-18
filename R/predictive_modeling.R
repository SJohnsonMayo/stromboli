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
