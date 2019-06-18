
getPermuteMatrixBlock <- function(permutations, strata) {
  strata <- factor(strata)
  strata <- factor(as.numeric(strata))
  res <- sapply(1:permutations, function(i) {
    strata1 <- strata
    levels(strata1) <- sample(levels(strata1))
    order(as.numeric(as.character(strata1)))
  })
  t(res)
}
# adonis2: permutate the covariate instead of data matrix - validated
# The variable of interest should be put at the last; Only p values of
# variable of interest will be reported Currently it can only be
# applied to subject-specific variable (nonvarying within subject) but
# covariates are allowed to vary Need to improve the speed
adonis2 <- function(formula, data = NULL, permutations = 999, method = "bray",
                    strata = NULL, block.perm = TRUE, contr.unordered = "contr.sum", contr.ordered = "contr.poly",
                    ...) {
  TOL <- 1e-07
  Terms <- terms(formula, data = data, keep.order = TRUE)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  qrhs <- qr(rhs)
  rhs <- rhs[, qrhs$pivot, drop = FALSE]
  rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
  grps <- grps[qrhs$pivot][1:qrhs$rank]
  u.grps <- unique(grps)
  nterms <- length(u.grps) - 1
  st <- ifelse(nterms == 1, 2, nterms)
  H.s <- lapply(st:(nterms + 1), function(j) {
    Xj <- rhs[, grps %in% u.grps[1:j]]
    qrX <- qr(Xj, tol = TOL)
    Q <- qr.Q(qrX)
    tcrossprod(Q[, 1:qrX$rank])
  })
  if (inherits(lhs, "dist")) {
    if (any(lhs < -TOL))
      stop("dissimilarities must be non-negative")
    dmat <- as.matrix(lhs^2)
  } else {
    dist.lhs <- as.matrix(vegdist(lhs, method = method, ...))
    dmat <- dist.lhs^2
  }
  n <- nrow(dmat)
  G <- -sweep(dmat, 1, rowMeans(dmat))/2
  SS.Exp.comb <- sapply(H.s, function(hat) sum(G * t(hat)))
  if (nterms == 1) {
    SS.Exp.each <- SS.Exp.comb
  } else {
    SS.Exp.each <- SS.Exp.comb[2] - SS.Exp.comb[1]
  }
  H.snterm <- H.s[[length(H.s)]]
  tIH.snterm <- t(diag(n) - H.snterm)
  # if (length(H.s) > 1) { H.s[[2]] <- H.s[[2]] - H.s[[1]] }
  SS.Res <- sum(G * tIH.snterm)
  df.Exp <- sum(grps == u.grps[length(u.grps)])
  df.Res <- n - qrhs$rank
  if (inherits(lhs, "dist")) {
    beta.sites <- qr.coef(qrhs, as.matrix(lhs))
    beta.spp <- NULL
  } else {
    beta.sites <- qr.coef(qrhs, dist.lhs)
    beta.spp <- qr.coef(qrhs, as.matrix(lhs))
  }
  colnames(beta.spp) <- colnames(lhs)
  colnames(beta.sites) <- rownames(lhs)
  F.Mod <- (SS.Exp.each/df.Exp)/(SS.Res/df.Res)
  f.test <- function(tH, G, df.Exp, df.Res, tIH.snterm) {
    (sum(G * tH)/df.Exp)/(sum(G * tIH.snterm)/df.Res)
  }
  rhs.1 <- rhs[, grps %in% u.grps[1:nterms], drop = F]
  rhs.2 <- rhs[, grps %in% u.grps[(nterms + 1)], drop = F]
  if (missing(strata))
    strata <- NULL
  if (block.perm == FALSE) {
    permute.ind <- vegan:::getPermuteMatrix(permutations, n, strata = strata)
  } else {
    if (is.null(strata))
      stop("Block permutation requires strata!\n")
    strata.u <- unique(strata)
    reorder.ind <- unlist(lapply(strata.u, function(x) which(strata ==
                                                               x)))
    expand.ind <- rep(1:length(strata.u), sapply(strata.u, function(x) sum(strata ==
                                                                             x)))
    rhs.2.u <- rhs.2[sapply(strata.u, function(x) which(strata == x)[1]),
                     , drop = F]
  }
  # To be checked and validated.
  f.perms <- sapply(1:permutations, function(i) {
    if (block.perm == FALSE) {
      rhs.2 <- rhs[permute.ind[i, ], grps %in% u.grps[(nterms + 1)],
                   drop = F]
    } else {
      rhs.2[reorder.ind, ] <- rhs.2.u[sample(nrow(rhs.2.u)), , drop = F][expand.ind,
                                                                         , drop = F]
    }
    Xj <- cbind(rhs.1, rhs.2)
    qrX <- qr(Xj, tol = TOL)
    Q <- qr.Q(qrX)
    tH.snterm <- t(tcrossprod(Q[, 1:qrX$rank]))
    tIH.snterm <- diag(n) - tH.snterm
    if (nterms > 1) {
      tH.snterm <- tH.snterm - t(H.s[[1]])
    }
    f.test(tH.snterm, G, df.Exp, df.Res, tIH.snterm)
  })
  f.perms <- round(f.perms, 12)
  F.Mod <- round(F.Mod, 12)
  SumsOfSqs = c(SS.Exp.each, SS.Res, SS.Exp.comb[length(SS.Exp.comb)] +
                  SS.Res)
  tab <- data.frame(Df = c(df.Exp, df.Res, n - 1), SumsOfSqs = SumsOfSqs,
                    MeanSqs = c(SS.Exp.each/df.Exp, SS.Res/df.Res, NA), F.Model = c(F.Mod,
                                                                                    NA, NA), R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)], P = c((sum(f.perms >=
                                                                                                                                                       F.Mod) + 1)/(permutations + 1), NA, NA))
  rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps[length(u.grps)]],
                     "Residuals", "Total")
  colnames(tab)[ncol(tab)] <- "Pr(>F)"
  attr(tab, "heading") <- "Terms added sequentially (first to last)\n"
  class(tab) <- c("anova", class(tab))
  out <- list(aov.tab = tab, call = match.call(), coefficients = beta.spp,
              coef.sites = beta.sites, f.perms = as.matrix(f.perms), model.matrix = rhs,
              terms = Terms)
  class(out) <- "adonis"
  out
}
######################################## This generates the matrix columns-wise From JnPaulson
generate_matrix <- function(x) {
  indptr = x$sample$matrix$indptr + 1
  indices = x$sample$matrix$indices + 1
  data = x$sample$matrix$data
  nr = length(x$observation$ids)
  counts = sapply(2:length(indptr), function(i) {
    x = rep(0, nr)
    seq = indptr[i - 1]:(indptr[i] - 1)
    x[indices[seq]] = data[seq]
    x
  })
  rownames(counts) = x$observation$ids
  colnames(counts) = x$sample$ids
  # I wish this next line wasn't necessary
  lapply(1:nrow(counts), function(i) {
    counts[i, ]
  })
}
generate_metadata <- function(x) {
  metadata = x$metadata
  metadata = lapply(1:length(x$ids), function(i) {
    id_metadata = lapply(metadata, function(j) {
      if (length(dim(j)) > 1) {
        as.vector(j[, i, drop = FALSE])
      } else {
        j[i]
      }
    })
    list(id = x$ids[i], metadata = id_metadata)
  })
  return(metadata)
}
namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(nm <- names(L)))
    nm <- snm
  if (any(nonames <- nm == ""))
    nm[nonames] <- snm[nonames]
  setNames(L, nm)
}
############################################ comm are otu counts, row: otus, column samples intersect.no: Pairwise
############################################ ratio calculated on pairs with at least 'intersect.no' common taxa;
############################################ Rev: 2016_12_12 add prevalence filter, which is not necessary Rev:
############################################ 2017_02_07
GMPR.old <- function(comm, intersect.no = 4, prev.filter = 0) {
  comm <- comm[rowMeans(comm != 0) >= prev.filter, ]
  ind.vec <- numeric(ncol(comm))
  res <- sapply(1:ncol(comm), function(i) {
    x <- comm[, i]
    pr <- sapply(1:ncol(comm), function(j) {
      y <- comm[, j]
      ind <- x != 0 & y != 0
      if (sum(ind) >= intersect.no) {
        res <- median(x[ind]/y[ind])
      } else {
        res <- NA
      }
    })
    if (sum(is.na(pr)) != 0)
      ind.vec[i] <<- 1
    exp(mean(log(pr[!is.na(pr)])))
  })
}
# New: 2017_02_07
GMPR <- function(comm, intersect.no = 4, ct.min = 2, trace = FALSE) {
  # Computes the GMPR size factor Args: comm: a matrix of counts, row -
  # features (OTUs, genes, etc) , column - sample intersect.no: the
  # minimum number of shared features between sample pair, where the
  # ratio is calculated ct.min: the minimum number of counts required to
  # calculate ratios ct.min = 5 has better results
  # Returns: a list that contains: gmpr： the GMPR size factors for all
  # samples; Samples with distinct sets of features will be output as NA.
  # nss: number of samples with significant sharing (> intersect.no)
  # including itself
  # mask counts < ct.min
  comm[comm < ct.min] <- 0
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0("S", 1:ncol(comm))
  }
  if (trace)
    cat("Begin GMPR size factor calculation ...\n")
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm), function(i) {
    if (i%%50 == 0) {
      if (trace)
        cat(i, "\n")
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x/comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm = TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  })
  if (sum(is.na(gmpr))) {
    warning(paste0("The following samples\n ", paste(colnames(comm)[is.na(gmpr)],
                                                     collapse = "\n"), "\ndo not share at least ", intersect.no,
                   " common taxa with the rest samples! ", "For these samples, their size factors are set to be NA! \n",
                   "You may consider removing these samples since they are potentially outliers or negative controls!\n",
                   "You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n"))
  }
  if (trace)
    cat("Completed!\n")
  if (trace)
    cat("Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n")
  attr(gmpr, "NSS") <- comm.no
  # Rev: 2017_09_07
  gmpr <- gmpr * median(colSums(comm))
  names(gmpr) <- colnames(comm)
  return(gmpr)
}
# New: 2018_02_24
dereplicate <- function(names) {
  names <- factor(names)
  ord <- order(order(names))
  tab <- table(names)
  suffix <- unlist(sapply(tab, function(x) {
    if (x == 1) {
      return("")
    } else {
      return(paste0("_", 1:x))
    }
  }))
  as.character(paste0(names, suffix[ord]))
}
is.na.null <- function(x) {
  if (is.null(x)) {
    return(TRUE)
  } else {
    if (is.na(x)[1]) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}
perform_demograph_analysis <- function(data.obj, grp.name) {
  obj <- summary(data.obj$meta.dat)
  write.csv(obj, "meta.data.summary.csv", quote = F)
  grp <- data.obj$meta.dat[, grp.name]
  res <- NULL
  pv.vec <- NULL
  if (is.factor(grp)) {
    for (var1 in setdiff(colnames(data.obj$meta.dat), grp.name)) {
      temp <- data.obj$meta.dat[, var1]
      res <- rbind(res, c("", var1, rep("", nlevels(grp) - 1)))
      res <- rbind(res, c("", levels(grp)))
      if (is.factor(temp)) {
        res <- rbind(res, cbind(levels(temp), table(temp, grp)))
        if (nlevels(temp) == 1) {
          pv <- NA
        } else {
          err <- try(pv <- formatC(fisher.test(table(temp, grp))$p.value,
                                   digit = 3))
          if (inherits(err, "try-error")) {
            pv <- NA
          }
        }
        res <- rbind(res, c("Fisher p", pv, rep("", nlevels(grp) -
                                                  1)))
      } else {
        res <- rbind(res, c("mean", aggregate(temp, by = list(grp),
                                              FUN = "mean")[, 2]))
        res <- rbind(res, c("sd", aggregate(temp, by = list(grp),
                                            FUN = "sd")[, 2]))
        err <- try(pv <- formatC(summary(aov(temp ~ grp))[[1]][1,
                                                               "Pr(>F)"], digit = 3))
        if (!inherits(err, "try-error")) {
          res <- rbind(res, c("ANOVA p", pv, rep("", nlevels(grp) -
                                                   1)))
        } else {
          res <- rbind(res, c("ANOVA p", "NA", rep("", nlevels(grp) -
                                                     1)))
          pv <- NA
        }
      }
      res <- rbind(res, rep("", nlevels(grp) + 1))
      pv.vec <- c(pv.vec, pv)
    }
  } else {
  }
  write.csv(res, "meta.data.by.grp.csv", row.names = F, quote = F)
  names(pv.vec) <- setdiff(colnames(data.obj$meta.dat), grp.name)
  return(pv.vec)
}
# New: 2018_02_25
estimate_richness_ <- function(OTU, measures = NULL) {
  # Modified over Phyloseq
  renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson",
                "Fisher", "Pielou")
  names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", "simpson",
                        "invsimpson", "fisher", "Pielou")
  if (is.null(measures)) {
    measures = as.character(renamevec)
  }
  if (any(measures %in% names(renamevec))) {
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in%
                                                            measures]
  }
  if (!any(measures %in% renamevec)) {
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  outlist = vector("list")
  estimRmeas = c("Chao1", "Observed", "ACE")
  if (any(estimRmeas %in% measures)) {
    outlist <- c(outlist, list(t(data.frame(vegan::estimateR(OTU)))))
  }
  if ("Shannon" %in% measures) {
    outlist <- c(outlist, list(shannon = vegan::diversity(OTU, index = "shannon")))
  }
  if ("Simpson" %in% measures) {
    outlist <- c(outlist, list(simpson = vegan::diversity(OTU, index = "simpson")))
  }
  if ("InvSimpson" %in% measures) {
    outlist <- c(outlist, list(invsimpson = vegan::diversity(OTU, index = "invsimpson")))
  }
  if ("Fisher" %in% measures) {
    fisher = tryCatch(vegan::fisher.alpha(OTU, se = TRUE), warning = function(w) {
      warning("estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
      suppressWarnings(vegan::fisher.alpha(OTU, se = TRUE)[, c("alpha",
                                                               "se")])
    })
    if (!is.null(dim(fisher))) {
      colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
      outlist <- c(outlist, list(fisher))
    } else {
      outlist <- c(outlist, Fisher = list(fisher))
    }
  }
  f <- function(ns) {
    sapply(ns, function(n) {
      i <- 1:n
      -sum(1/i * log(1/i))
    })
  }
  cat(names(outlist))
  if ("Pielou" %in% measures) {
    outlist <- c(outlist, list(Pielou = vegan::diversity(OTU, index = "shannon")/f(rowSums(OTU !=
                                                                                             0))))
  }
  out = do.call("cbind", outlist)
  namechange = intersect(colnames(out), names(renamevec))
  colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
  colkeep = sapply(paste0("(se\\.){0,}", measures), grep, colnames(out),
                   ignore.case = TRUE)
  out = out[, sort(unique(unlist(colkeep))), drop = FALSE]
  out <- as.data.frame(out)
  return(out)
}
# Rev: 2016_09_10, Implement p value based omnibus test
PermanovaG2 <- function(formula, dat = NULL, ...) {
  save.seed <- get(".Random.seed", .GlobalEnv)
  lhs <- formula[[2]]
  lhs <- eval(lhs, dat, parent.frame())
  rhs <- as.character(formula)[3]
  p.perms <- list()
  p.obs <- list()
  for (i in 1:(dim(lhs)[3])) {
    assign(".Random.seed", save.seed, .GlobalEnv)
    Y <- as.dist(lhs[, , i])
    formula2 <- as.formula(paste("Y", "~", rhs))
    obj <- adonis(formula2, dat, ...)
    perm.mat <- obj$f.perms
    p.perms[[i]] <- 1 - (apply(perm.mat, 2, rank) - 1)/nrow(perm.mat)
    p.obs[[i]] <- obj$aov.tab[1:ncol(perm.mat), "Pr(>F)"]
  }
  omni.pv <- NULL
  indiv.pv <- NULL
  for (j in 1:ncol(perm.mat)) {
    p.perms.j <- sapply(p.perms, function(x) x[, j])
    p.obj.j <- sapply(p.obs, function(x) x[j])
    omni.pv <- c(omni.pv, mean(c(rowMins(p.perms.j) <= min(p.obj.j),
                                 1)))
    indiv.pv <- rbind(indiv.pv, p.obj.j)
  }
  colnames(indiv.pv) <- paste0("D", 1:ncol(indiv.pv), ".p.value")
  rownames(indiv.pv) <- 1:nrow(indiv.pv)
  aov.tab <- data.frame(indiv.pv, omni.p.value = omni.pv)
  rownames(aov.tab) <- rownames(obj$aov.tab)[1:ncol(perm.mat)]
  list(aov.tab = aov.tab)
}
distance_compare_test <- function(dist.mat, ind1, ind2, ind3, ID2 = NULL,
                                  ID3 = NULL, alternative = "greater", nperm = 999) {
  ind23 <- c(ind2, ind3)
  if (!is.null(ID2) & !is.null(ID3)) {
    ID23 <- c(ID2, ID3)
    ID2u <- unique(ID2)
    ID3u <- unique(ID3)
    ID23u <- c(ID2u, ID3u)
    n2u <- length(ID2u)
    n3u <- length(ID3u)
  }
  n2 <- length(ind2)
  n3 <- length(ind3)
  dist12 <- dist.mat[ind1, ind2]
  dist13 <- dist.mat[ind1, ind3]
  stat.obs <- mean(dist12) - mean(dist13)
  stat.perm <- sapply(1:nperm, function(i) {
    if (is.null(ID2) & is.null(ID3)) {
      ind23.p <- sample(ind23)
      ind2.p <- ind23.p[1:n2]
      ind3.p <- ind23.p[(n2 + 1):(n2 + n3)]
    } else {
      ID23u.p <- sample(ID23u)
      ID2.p <- ID23u.p[1:n2u]
      ID3.p <- ID23u.p[(n2u + 1):(n2u + n3u)]
      ind2.p <- ind23[ID23 %in% ID2.p]
      ind3.p <- ind23[ID23 %in% ID3.p]
    }
    dist12.p <- dist.mat[ind1, ind2.p]
    dist13.p <- dist.mat[ind1, ind3.p]
    mean(dist12.p) - mean(dist13.p)
  })
  if (alternative == "greater") {
    pv <- mean(c(stat.perm >= stat.obs, TRUE))
  } else {
    pv <- mean(c(stat.perm <= stat.obs, TRUE))
  }
  pv
}
twopart.test <- function(x1, x2, zero.p = 0.2) {
  n1 <- length(x1)
  n2 <- length(x2)
  p1 <- mean(x1 != 0)
  p2 <- mean(x2 != 0)
  m1 <- sum(x1 != 0)
  m2 <- sum(x2 != 0)
  p12 <- (m1 + m2)/(n1 + n2)
  q12 <- 1 - p12
  if (q12 >= zero.p) {
    Z <- (abs(p1 - p2) - (1/(2 * n1) + 1/(2 * n2)))/sqrt(p12 * q12 *
                                                           (1/n1 + 1/n2))
    x1 <- x1[x1 != 0]
    x2 <- x2[x2 != 0]
    R1 <- sum(rank(c(x1, x2))[1:length(x1)])
    ti <- as.vector(table(c(x1, x2)))
    W <- (abs(R1 - m1 * (m1 + m2 + 1)/2) - 1/2)/sqrt((m1 * m2/12) *
                                                       (m1 + m2 + 1 - sum(ti * (ti^2 - 1))/(m1 + m2)/(m1 + m2 - 1)))
    X2 <- Z^2 + W^2
    res <- list()
    res$stat <- X2
    res$p.value <- 1 - pchisq(X2, 2)
    res$Z <- Z
    res$W <- W
    res$test <- "TwoPart"
  } else {
    res <- wilcox.test(x1, x2)
    res$test <- "Wilcox"
  }
  res
}
getPermuteMatrix <- function(perm, N, strata = NULL) {
  if (length(perm) == 1) {
    perm <- how(nperm = perm)
  }
  if (!missing(strata) && !is.null(strata)) {
    if (inherits(perm, "how") && is.null(getBlocks(perm)))
      setBlocks(perm) <- strata
  }
  if (inherits(perm, "how"))
    perm <- shuffleSet(N, control = perm)
  if (is.null(attr(perm, "control")))
    attr(perm, "control") <- structure(list(within = list(type = "supplied matrix"),
                                            nperm = nrow(perm)), class = "how")
  perm
}
perm_fdr_adj <- function(F0, Fp) {
  ord <- order(F0, decreasing = T)
  F0 <- F0[ord]
  perm.no <- ncol(Fp)
  Fp <- as.vector(Fp)
  Fp <- Fp[!is.na(Fp)]
  Fp <- sort(c(Fp, F0), decreasing = F)
  n <- length(Fp)
  m <- length(F0)
  FPN <- (n + 1) - match(F0, Fp) - 1:m
  p.adj.fdr <- FPN/perm.no/(1:m)
  # p.adj.fdr <- sapply(F0, function(x) sum(Fp >= x, na.rm=TRUE) /
  # perm.no)/(1:length(F0))
  p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
}
perm_fwer_adj <- function(F0, Fp) {
  ord <- order(F0, decreasing = T)
  F0 <- F0[ord]
  col.max <- colMaxs(Fp, na.rm = TRUE)
  p.adj.fwer <- sapply(F0, function(x) mean(col.max >= x))[order(ord)]
}
# New: 2018_06_20 Westfall young
perm_fwer_adj2 <- function(F0, Fp) {
  ord <- order(F0, decreasing = T)
  m <- length(F0)
  F0 <- F0[ord]
  Fp <- Fp[ord, , drop = FALSE]
  col.max <- Fp[m, ]
  p.adj.fwer <- sapply(m:1, function(i) {
    x <- F0[i]
    y <- Fp[i, ]
    col.max <<- ifelse(y > col.max, y, col.max)
    mean(col.max >= x)
  })
  p.adj.fwer <- rev(p.adj.fwer)
  p.adj.fwer <- pmin(1, rev(cummin(rev(p.adj.fwer))))[order(ord)]
}
# Need to be further comprehensively tested a. Speed up, finished!  b.
# Permutation method （response, covariate or residual permutation),
# residual will be default!  c. Revise - add permutation stratified by
# subject, finished!  Rev: 2017_02_02 permutation-based FDR control
# Rev: 2017_02_13 Add LMM-based permutation for block.perm=TRUE (type I
# error is controled, but power study hasn't been comprehensively
# studied) Rev: 2017_02_24 allow 'adj.name' to contain multiple
# covariates Still need to address NA's, currently simply remove NA's
# Rev: 2018_03_18 add winsorization support for permutation test - Only
# one side winsorization Rev: 2018_03_18 add effect size (partial R2)
# and Direction Rev: 2018_06_20 add Westfall young FWER adjustment
permute_differential_analysis <- function(meta.dat, prop, grp.name, adj.name = NULL,
                                          strata = NULL, block.perm = FALSE, sqrt.trans = TRUE, resid.perm = TRUE,
                                          perm.no = 999, winsor = FALSE, winsor.qt = 0.97, trace = FALSE) {
  # Square root transformation User should take care of the
  # normalization, transformation and addressing outliers
  opt.old <- options("contrasts")
  options(contrasts = c("contr.treatment", "contr.poly"))
  meta.dat <- droplevels(meta.dat)
  if (sqrt.trans) {
    Y <- sqrt(prop)
  } else {
    Y <- prop
  }
  row.names <- rownames(Y)
  if (winsor == TRUE) {
    # Only one side winsorization - be careful for different type of
    # transformation
    Y <- t(apply(Y, 1, function(x) {
      cutoff <- quantile(x, winsor.qt)
      x[x >= cutoff] <- cutoff
      x
    }))
  }
  if (!is.null(strata)) {
    strata <- factor(strata)
  }
  # Prepare model matrix
  n <- ncol(prop)
  I <- diag(n)
  if (is.null(adj.name)) {
    M0 <- model.matrix(~1, meta.dat)
  } else {
    df0 <- meta.dat[, c(adj.name), drop = F]
    M0 <- model.matrix(~., df0)
  }
  # P0 <- M0 %*% solve(t(M0) %*% M0) %*% t(M0)
  df1 <- meta.dat[, c(adj.name, grp.name), drop = F]
  M1 <- model.matrix(~., df1)
  # Make sure it is contr.tretament
  M2 <- model.matrix(~., meta.dat[, c(grp.name), drop = F])
  # Remove the effect of confounder
  M2 <- as.matrix(resid(lm(M2[, -1] ~ M0 - 1)))
  # QR decompostion
  qrX0 <- qr(M0, tol = 1e-07)
  Q0 <- qr.Q(qrX0)
  Q0 <- Q0[, 1:qrX0$rank, drop = FALSE]
  qrX1 <- qr(M1, tol = 1e-07)
  Q1 <- qr.Q(qrX1)
  Q1 <- Q1[, 1:qrX1$rank, drop = FALSE]
  # Got residual
  if (resid.perm) {
    if (!block.perm) {
      # Permute the residual
      if (is.null(adj.name)) {
        Y <- t(resid(lm(as.formula(paste("t(Y) ~ 1")), meta.dat)))
      } else {
        Y <- t(resid(lm(as.formula(paste("t(Y) ~ ", paste(adj.name,
                                                          collapse = "+"))), meta.dat)))
      }
    } else {
      if (is.null(strata)) {
        stop("Block permutation requires strata!\n")
      } else {
        Y.r <- matrix(NA, nrow(Y), nlevels(strata))
        Y.e <- Y
        if (trace)
          cat("Fitting linear mixed effects model ...\n")
        for (j in 1:nrow(Y)) {
          # Linear mixed effects model
          yy <- Y[j, ]
          meta.dat$yy <- yy
          meta.dat$strata <- strata
          if (is.null(adj.name)) {
            # obj <- lme(as.formula(paste('yy ~ 1')), random =~ 1 | strata,
            # data=meta.dat, method='ML')
            obj <- lmer(as.formula(paste("yy ~ 1 + (1 | strata)")),
                        data = meta.dat, REML = FALSE)
            # The order is the same as the levels Y.r[j, ] <- random.effects(obj)[,
            # 1]
            Y.r[j, ] <- random.effects(obj)[[1]][, 1]
            Y.e[j, ] <- resid(obj)
          } else {
            # obj <- lme(as.formula(paste('yy ~ ', paste(adj.name, collapse='+'))),
            # random =~ 1 | strata, data=meta.dat, method='ML')
            obj <- lmer(as.formula(paste("yy ~ ", paste(adj.name,
                                                        collapse = "+"), " + (1 | strata)")), data = meta.dat,
                        REML = FALSE)
            # Y.r[j, ] <- random.effects(obj)[, 1]
            Y.r[j, ] <- random.effects(obj)[[1]][, 1]
            Y.e[j, ] <- resid(obj)
          }
        }
        Y <- Y.r[, as.numeric(strata)] + Y.e
        # Y <- Y - rowMeans(Y) Y <- Y.e
      }
    }
  }
  TSS <- rowSums(Y^2)
  MSS1 <- rowSums((Y %*% Q1)^2)
  MSS0 <- rowSums((Y %*% Q0)^2)  # Not necessary, it's zero
  F0 <- (MSS1 - MSS0)/(TSS - MSS1)
  R2 <- (MSS1 - MSS0)/(TSS - MSS0)
  coefficients <- t(solve(t(M2) %*% M2) %*% t(M2) %*% t(Y))
  # P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1) F0 <- diag(Y %*% (P1 - P0)
  # %*% t(Y)) / diag(Y %*% (I - P1) %*% t(Y)) df3 <- df1
  if (block.perm == FALSE) {
    perm.ind <- vegan:::getPermuteMatrix(perm.no, n, strata = strata)
    perm.no <- nrow(perm.ind)
  }
  if (trace)
    cat("Permutation test ....\n")
  Fp <- sapply(1:perm.no, function(i) {
    if (i%%100 == 0)
      if (trace)
        cat(".")
    if (block.perm == FALSE) {
      Yp <- Y[, perm.ind[i, ]]
    } else {
      # Double permutation
      strata.p <- factor(strata, levels = sample(levels(strata)))
      Yp <- Y.r[, as.numeric(strata.p)] + Y.e[, sample(ncol(Y))]
      # Yp <- Y.e[, sample(ncol(Y))] Yp <- Yp - rowMeans(Yp)
    }
    # df3[, grp.name] <- sample(df1[, grp.name]) M1 <- model.matrix( ~.,
    # df3) P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
    MSS1p <- rowSums((Yp %*% Q1)^2)
    MSS0p <- rowSums((Yp %*% Q0)^2)
    if (block.perm == FALSE) {
      TSSp <- TSS
    } else {
      TSSp <- rowSums(Yp^2)
    }
    (MSS1p - MSS0p)/(TSSp - MSS1p)
  })
  if (mean(is.na(F0)) >= 0.1) {
    warning("More than 10% observed F stats have NA! Please check! \n")
  }
  if (mean(is.na(Fp)) >= 0.1) {
    warning("More than 10% permuted F stats have NA! Please check! \n")
  }
  na.ind <- is.na(F0)
  F0 <- F0[!na.ind]
  Fp <- Fp[!na.ind, ]
  p.raw <- cbind(Fp >= F0, 1)
  p.raw <- rowMeans(p.raw)
  # p.raw[is.na(p.raw)] <- 1
  p.adj.fdr <- perm_fdr_adj(F0, Fp)
  p.adj.fwer <- perm_fwer_adj2(F0, Fp)
  # Pad back the NA values
  pad <- function(vec, ind) {
    vec0 <- numeric(length(ind))
    vec0[!ind] <- vec
    vec0[ind] <- NA
    vec0
  }
  options(contrasts = opt.old$contrast)
  F0 <- pad(F0, na.ind)
  p.raw <- pad(p.raw, na.ind)
  p.adj.fdr <- pad(p.adj.fdr, na.ind)
  p.adj.fwer <- pad(p.adj.fwer, na.ind)
  names(F0) <- names(p.raw) <- names(p.adj.fdr) <- names(p.adj.fwer) <- row.names
  return(list(F.stat = F0, R2 = R2, coefficients = coefficients, p.raw = p.raw,
              p.adj.fdr = p.adj.fdr, p.adj.fwer = p.adj.fwer))
  # return(p.raw)
}
#
# New: 2017_08_17 Add a new variant of PERMANOVA with matrix
# decomposition Partially validated
permanova2 <- function(meta.dat, D, grp.name, adj.name = NULL, strata = NULL,
                       block.perm = FALSE, resid.perm = TRUE, perm.no = 999, eig = c("All",
                                                                                     "Positive", "Top"), var.exp = 0.9) {
  eig <- match.arg(eig)
  # Square root transformation User should take care of the
  # normalization, transformation and addressing outliers
  D <- -as.matrix(D)^2/2
  n <- nrow(D)
  G <- mean(D) + D - rowMeans(D) - matrix(rep(1, n), ncol = 1) %*% colMeans(D)
  eig.obj <- eigen(G)
  if (eig == "All") {
    Y <- eig.obj$vectors
    lambda <- eig.obj$values
    ind <- abs(lambda) > 1e-06
    lambda <- lambda[ind]
    Y <- Y[, ind, drop = FALSE]
  }
  if (eig == "Positive") {
    lambda <- eig.obj$values
    ind <- lambda > 1e-06
    lambda <- lambda[ind]
    Y <- Y[, ind, drop = FALSE]
  }
  if (eig == "Top") {
    Y <- eig.obj$vectors
    lambda <- eig.obj$values
    ind <- lambda > 1e-06
    lambda <- lambda[ind]
    Y <- Y[, ind, drop = FALSE]
    ind <- 1:(which(cumsum(lambda/sum(lambda)) > var.exp)[1])
    lambda <- lambda[ind]
    Y <- Y[, ind, drop = FALSE]
  }
  Y <- t(Y)
  if (!is.null(strata)) {
    strata <- factor(strata)
  }
  # Prepare model matrix
  n <- ncol(Y)
  I <- diag(n)
  if (is.null(adj.name)) {
    M0 <- model.matrix(~1, meta.dat)
  } else {
    df0 <- meta.dat[, c(adj.name), drop = F]
    M0 <- model.matrix(~., df0)
  }
  # P0 <- M0 %*% solve(t(M0) %*% M0) %*% t(M0)
  df1 <- meta.dat[, c(adj.name, grp.name), drop = F]
  M1 <- model.matrix(~., df1)
  # QR decompostion
  qrX0 <- qr(M0, tol = 1e-07)
  Q0 <- qr.Q(qrX0)
  Q0 <- Q0[, 1:qrX0$rank, drop = FALSE]
  qrX1 <- qr(M1, tol = 1e-07)
  Q1 <- qr.Q(qrX1)
  Q1 <- Q1[, 1:qrX1$rank, drop = FALSE]
  # Got residual
  if (resid.perm) {
    if (!block.perm) {
      # Permute the residual
      if (is.null(adj.name)) {
        Y <- t(resid(lm(as.formula(paste("t(Y) ~ 1")), meta.dat)))
      } else {
        Y <- t(resid(lm(as.formula(paste("t(Y) ~ ", paste(adj.name,
                                                          collapse = "+"))), meta.dat)))
      }
    } else {
      if (is.null(strata)) {
        stop("Block permutation requires strata!\n")
      } else {
        Y.r <- matrix(NA, nrow(Y), nlevels(strata))
        Y.e <- Y
        # cat('Fitting linear mixed effects model ...\n')
        for (j in 1:nrow(Y)) {
          # Linear mixed effects model
          yy <- Y[j, ]
          meta.dat$yy <- yy
          meta.dat$strata <- strata
          if (is.null(adj.name)) {
            obj <- lmer(as.formula(paste("yy ~ 1 + (1 | strata)")),
                        data = meta.dat)
            # The order is the same as the levels
            Y.r[j, ] <- random.effects(obj)[[1]][, 1]
            Y.e[j, ] <- resid(obj)
          } else {
            obj <- lmer(as.formula(paste("yy ~ ", paste(adj.name,
                                                        collapse = "+"), "+ (1 | strata)")), data = meta.dat)
            Y.r[j, ] <- random.effects(obj)[[1]][, 1]
            Y.e[j, ] <- resid(obj)
          }
        }
        Y <- Y.r[, as.numeric(strata)] + Y.e
        # Y <- Y - rowMeans(Y) Y <- Y.e
      }
    }
  }
  TSS <- rowSums(Y^2)
  MSS1 <- rowSums((Y %*% Q1)^2)
  MSS0 <- rowSums((Y %*% Q0)^2)  # Not necessary, it's zero
  F0 <- sum(lambda * (MSS1 - MSS0))/sum(lambda * (TSS - MSS1))
  # P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1) F0 <- diag(Y %*% (P1 - P0)
  # %*% t(Y)) / diag(Y %*% (I - P1) %*% t(Y)) df3 <- df1
  if (block.perm == FALSE) {
    perm.ind <- vegan:::getPermuteMatrix(perm.no, n, strata = strata)
    perm.no <- nrow(perm.ind)
  }
  # cat('Permutation test ....\n')
  Fp <- sapply(1:perm.no, function(i) {
    # if (i %% 100 == 0) cat('.')
    if (block.perm == FALSE) {
      Yp <- Y[, perm.ind[i, ]]
    } else {
      # Double permutation
      strata.p <- factor(strata, levels = sample(levels(strata)))
      Yp <- Y.r[, as.numeric(strata.p)] + Y.e[, sample(ncol(Y))]
      # Yp <- Y.e[, sample(ncol(Y))] Yp <- Yp - rowMeans(Yp)
    }
    # df3[, grp.name] <- sample(df1[, grp.name]) M1 <- model.matrix( ~.,
    # df3) P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
    MSS1p <- rowSums((Yp %*% Q1)^2)
    MSS0p <- rowSums((Yp %*% Q0)^2)
    if (block.perm == FALSE) {
      TSSp <- TSS
    } else {
      TSSp <- rowSums(Yp^2)
    }
    sum(lambda * (MSS1p - MSS0p))/sum(lambda * (TSSp - MSS1p))
  })
  p.value <- mean(c(Fp >= F0, 1))
  return(list(f0 = F0, f.perms = Fp, p.value = p.value))
  # return(p.raw)
}
# This function for nonparametric/permutaiton method Rev: 2017_02_16
# Add normalization method; Add transformation; Remove rarefaction
# (only output warnings); Rev: 2017_10_30 Support filtering based on
# coefficient of variation Rev: 2018_01_30 Add ct.min and handle NA
# size factor Rev: 2018_03_15 Add winsorization support for permutatin
# test - default is FALSE Rev: 2018_10_03 TSS size factor calculation
# is based on specific taxa level
perform_differential_analysis <- function(data.obj, grp.name, adj.name = NULL,
                                          subject = NULL, taxa.levels = c("Phylum", "Order", "Class", "Family",
                                                                          "Genus", "Species"), method = "perm", block.perm = FALSE, perm.no = 999,
                                          norm = "GMPR", norm.level = "Species", intersect.no = 4, ct.min = 2,
                                          transform = "sqrt", winsor = FALSE, winsor.qt = 0.97, prev = 0.1, minp = 0.002,
                                          medianp = NULL, cv = NULL, mt.method = "fdr", cutoff = 0.15, ann = "",
                                          seed = 123, ...) {
  # To be completed
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  ind <- !is.na(grp)
  data.obj <- subset_data(data.obj, ind)
  grp <- grp[ind]
  df <- df[ind, ]
  if ("Species" %in% taxa.levels & !("Species" %in% names(data.obj$abund.list))) {
    data.obj$abund.list[["Species"]] <- data.obj$otu.tab
    rownames(data.obj$abund.list[["Species"]]) <- paste0("OTU", rownames(data.obj$otu.tab),
                                                         ":", data.obj$otu.name[, "Phylum"], ";", data.obj$otu.name[,
                                                                                                                    "Genus"])
  }
  # Test for sequence-depth confounding
  dep <- colSums(data.obj$otu.tab)
  diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, "Pr(>F)"]
  if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
    warning(paste0("\nSignificant sequencing depth confounding with the variable of interest!\n",
                   "For nonparametric test/permutaiton test, there may be potentially many false postives (for those less prevalent taxa)!\n",
                   "Consider performing rarefaction first! (However, rarefaction will not completly solve the problem.)\n",
                   "May also try count-based models, which might have better false postive controls!\n"))
    # cat('For parametric test with sequence depth adjustment (DESeq2),
    # please be cautious about the results!\n There may be potential
    # residual sequence depth confounding!\n') potential revising
  }
  # Calculate size.factor
  if (norm == "Precalculated") {
    size.factor <- data.obj$size.factor
  }
  if (norm == "GMPR") {
    if (norm.level %in% c("OTU", "Species")) {
      tab <- data.obj$otu.tab
    } else {
      tab <- data.obj$abund.list[[norm.level]]
    }
    size.factor <- GMPR(tab, intersect.no, ct.min)
    # Rev: 2018_01_30
    ind <- !is.na(size.factor)
    data.obj <- subset_data(data.obj, ind)
    grp <- grp[ind]
    df <- df[ind, ]
    size.factor <- size.factor[ind]
  }
  # Method-dependent processing
  if (is.null(method)) {
    if (nlevels(grp) == 2) {
      method <- "wilcox"
    } else {
      method <- "kruskal"
    }
  }
  if (method == "wilcox.pair") {
    if (nlevels(grp) != 2)
      stop("Wilcox test requires two groups!\n")
    if (is.null(subject))
      stop("Paired wilcox needs subject information!\n")
    subject <- factor(df[, subject])
    ind1 <- ind2 <- NULL
    for (sub in levels(subject)) {
      temp1 <- which(as.numeric(grp) == 1 & subject == sub)
      temp2 <- which(as.numeric(grp) == 2 & subject == sub)
      if (length(temp1) != 0 & length(temp2) != 0) {
        ind1 <- c(ind1, temp1[1])
        ind2 <- c(ind2, temp2[1])
      }
    }
    # if (length(ind1) != length(ind2)) warning('Some subjects are not
    # paired!\n')
  }
  if (method == "perm.pair") {
    if (is.null(subject))
      stop("Paired permutation test needs subject information!\n")
    subject <- factor(df[, subject])
  }
  if (method == "perm") {
    if (!is.null(subject)) {
      subject <- factor(df[, subject])
    }
  }
  ann <- paste0(ann, ".", method)
  if (is.factor(grp)) {
    pv.list <- qv.list <- qv.perm.list <- R2.list <- coef.list <- fc.list <- pc.list <- m.list <- nzm.list <- prv.list <- list()
    res.final <- NULL
    for (LOI in taxa.levels) {
      cat(LOI, "\n")
      ct <- data.obj$abund.list[[LOI]]
      # Rev: 2018_10_03 Move here
      if (norm == "TSS") {
        size.factor <- colSums(data.obj$abund.list[[LOI]])
      }
      # Filtering
      prop0 <- t(t(ct)/colSums(ct))
      if (!is.null(prev)) {
        prop0 <- prop0[rowSums(prop0 != 0) > prev * ncol(prop0),
                       , drop = FALSE]
        ct <- ct[rownames(prop0), , drop = FALSE]
      }
      if (!is.null(minp)) {
        prop0 <- prop0[rowMaxs(prop0) > minp, , drop = FALSE]
        ct <- ct[rownames(prop0), , drop = FALSE]
      }
      if (!is.null(medianp)) {
        nz.mean <- apply(prop0, 1, function(x) median(x[x != 0]))
        prop0 <- prop0[nz.mean > medianp, , drop = FALSE]
        ct <- ct[rownames(prop0), , drop = FALSE]
      }
      if (!is.null(cv)) {
        prop0 <- prop0[rowSds(prop0)/rowMeans(prop0) > cv, , drop = FALSE]
        ct <- ct[rownames(prop0), , drop = FALSE]
      }
      # Normalization
      prop <- t(t(ct)/size.factor)
      # Transformation - Others are possible/may be explored in the future
      if (transform == "sqrt") {
        prop <- sqrt(prop)
      }
      if (method == "perm") {
        set.seed(seed)
        perm.obj <- permute_differential_analysis(df, prop, grp.name,
                                                  adj.name, strata = subject, block.perm = block.perm,
                                                  sqrt.trans = FALSE, perm.no = perm.no, winsor = winsor,
                                                  winsor.qt = winsor.qt)
        pv.de2 <- perm.obj$p.raw
        names(pv.de2) <- rownames(prop)
      }
      # For legacy use
      if (method == "perm.pair") {
        set.seed(seed)
        perm.obj <- permute_differential_analysis(df, prop, grp.name,
                                                  adj.name, strata = subject, block.perm = FALSE, sqrt.trans = FALSE,
                                                  perm.no = perm.no, winsor = winsor, winsor.qt = winsor.qt)
        pv.de2 <- perm.obj$p.raw
        names(pv.de2) <- rownames(prop)
      }
      pv.vec <- m.vec <- nzm.vec <- prv.vec <- fc.vec <- pc.vec <- NULL
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        taxon.abund0 <- prop0[taxon, ]
        pv <- fc <- pc <- m <- nzm <- prv <- NULL
        if (method == "wilcox") {
          if (nlevels(grp) != 2)
            stop("Wilcox test requires two groups!\n")
          pv <- wilcox.test(taxon.abund ~ grp)$p.value
        }
        if (method == "wilcox.pair") {
          pv <- wilcox.test(taxon.abund[ind1], taxon.abund[ind2],
                            paired = T)$p.value
        }
        if (method == "twopart") {
          if (nlevels(grp) != 2)
            stop("Two part test requires two groups!\n")
          grp1 <- taxon.abund[as.numeric(grp) == 1]
          grp2 <- taxon.abund[as.numeric(grp) == 2]
          pv <- twopart.test(grp1, grp2)$p.value
        }
        if (method == "kruskal") {
          if (nlevels(grp) <= 2)
            warning("Kruskal-wallis test requires three or more groups!\n")
          pv <- kruskal.test(taxon.abund ~ grp)$p.value
        }
        if (method == "perm") {
          pv <- pv.de2[taxon]
        }
        if (method == "perm.pair") {
          pv <- pv.de2[taxon]
        }
        m <- tapply(taxon.abund0, grp, function(x) mean(x))
        nzm <- tapply(taxon.abund0, grp, function(x) mean(x[x !=
                                                              0]))
        prv <- tapply(taxon.abund0, grp, function(x) sum(x != 0))
        # Rev: 2017_02_21 fc change baseline grp
        if (nlevels(grp) == 2) {
          grp.no <- table(grp)
          fc <- log2(m[2]/m[1])
          pc <- prv[2]/grp.no[2]/prv[1] * grp.no[1]
        } else {
          pc <- fc <- NA
        }
        pv.vec <- rbind(pv.vec, pv)
        fc.vec <- rbind(fc.vec, fc)
        m.vec <- rbind(m.vec, m)
        nzm.vec <- rbind(nzm.vec, nzm)
        pc.vec <- rbind(pc.vec, pc)
        prv.vec <- rbind(prv.vec, prv/table(grp))
      }
      temp <- p.adjust(pv.vec[, 1], "fdr")
      qv.vec <- matrix(temp, ncol = 1)
      rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(pc.vec) <- rownames(m.vec) <- rownames(nzm.vec) <- rownames(prv.vec) <- rownames(prop)
      colnames(pv.vec) <- "Pvalue"
      colnames(qv.vec) <- "Qvalue"
      colnames(fc.vec) <- "log2(AbundRatio)"
      colnames(pc.vec) <- "PrevRatio"
      colnames(m.vec) <- paste(levels(grp), "Mean")
      colnames(nzm.vec) <- paste(levels(grp), "nzMean")
      colnames(prv.vec) <- paste(levels(grp), "preval")
      pv.list[[LOI]] <- pv.vec
      qv.list[[LOI]] <- qv.vec
      fc.list[[LOI]] <- fc.vec
      pc.list[[LOI]] <- pc.vec
      m.list[[LOI]] <- m.vec
      nzm.list[[LOI]] <- nzm.vec
      prv.list[[LOI]] <- prv.vec
      if (method %in% c("perm", "perm.pair")) {
        qv.perm.list[[LOI]] <- t(t(perm.obj$p.adj.fdr))
        R2.list[[LOI]] <- t(t(perm.obj$R2))
        coef.list[[LOI]] <- perm.obj$coefficients
      }
      if (method %in% c("perm", "perm.pair")) {
        coef.mat <- perm.obj$coefficients
        colnames(coef.mat) <- paste0("coef_", colnames(coef.mat))
        res <- cbind(pv.vec, qv.vec, Qvalue.perm = perm.obj$p.adj.fdr,
                     R2 = perm.obj$R2, coef.mat, m.vec, nzm.vec, fc.vec, prv.vec,
                     pc.vec)
        rownames(res) <- rownames(prop)
      } else {
        res <- cbind(pv.vec, qv.vec, m.vec, nzm.vec, fc.vec, prv.vec,
                     pc.vec)
        rownames(res) <- rownames(prop)
      }
      write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_",
                            LOI, "_", ann, ".csv"))
      if (mt.method == "fdr") {
        res.final <- rbind(res.final, res[res[, "Qvalue"] <= cutoff,
                                          , drop = F])
      }
      if (mt.method == "raw") {
        res.final <- rbind(res.final, res[res[, "Pvalue"] <= cutoff,
                                          , drop = F])
      }
      if (mt.method == "fdr.perm") {
        res.final <- rbind(res.final, res[res[, "Qvalue.perm"] <=
                                            cutoff, , drop = F])
      }
    }
    if (!is.null(res.final)) {
      colnames(res.final) <- colnames(res)
      write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_",
                                  mt.method, "_", cutoff, "_", ann, ".csv"))
    }
    return(list(pv.list = pv.list, fc.list = fc.list, pc.list = pc.list,
                qv.list = qv.list, qv.perm.list = qv.perm.list, R2.list = R2.list,
                coef.list = coef.list, m.list = m.list))
  } else {
    if (is.null(method)) {
      method <- "Spearman"
    }
    # Continuous case - currently only has DESeq2
    pv.list <- qv.list <- fc.list <- m.list <- R2.list <- coef.list <- qv.perm.list <- list()
    res.final <- NULL
    for (LOI in taxa.levels) {
      cat(LOI, "\n")
      ct <- data.obj$abund.list[[LOI]]
      # Rev: 2018_12_25
      if (norm == "TSS") {
        size.factor <- colSums(data.obj$abund.list[[LOI]])
      }
      # Filtering
      prop0 <- t(t(ct)/colSums(ct))
      if (!is.null(prev)) {
        prop0 <- prop0[rowSums(prop0 != 0) > prev * ncol(prop0),
                       , drop = FALSE]
        ct <- ct[rownames(prop0), , drop = FALSE]
      }
      if (!is.null(minp)) {
        prop0 <- prop0[rowMaxs(prop0) > minp, , drop = FALSE]
        ct <- ct[rownames(prop0), , drop = FALSE]
      }
      if (!is.null(medianp)) {
        nz.mean <- apply(prop0, 1, function(x) median(x[x != 0]))
        prop0 <- prop0[nz.mean > medianp, , drop = FALSE]
        ct <- ct[rownames(prop0), , drop = FALSE]
      }
      if (!is.null(cv)) {
        prop0 <- prop0[rowSds(prop0)/rowMeans(prop0) > cv, , drop = FALSE]
        ct <- ct[rownames(prop0), , drop = FALSE]
      }
      # Normalization
      prop <- t(t(ct)/size.factor)
      # Transformation - Others are possible/may be explored in the future
      if (transform == "sqrt") {
        prop <- sqrt(prop)
      }
      if (method == "perm") {
        set.seed(seed)
        perm.obj <- permute_differential_analysis(df, prop, grp.name,
                                                  adj.name, strata = subject, block.perm = block.perm,
                                                  sqrt.trans = FALSE, perm.no = perm.no, winsor = winsor,
                                                  winsor.qt = winsor.qt)
        pv.de2 <- perm.obj$p.raw
        names(pv.de2) <- rownames(prop)
        # Place holder
        fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i,
                                                            ], df[, grp.name], method = "spearman"))
      }
      if (method == "Spearman") {
        if (!is.null(adj.name)) {
          stop("Spearman test can't adjust covariates!")
        }
        pv.de2 <- apply(prop, 1, function(x) {
          cor.test(x, grp, method = "spearman")$p.value
        })
        names(pv.de2) <- rownames(prop)
        # Place holder
        fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i,
                                                            ], df[, grp.name], method = "spearman"))
      }
      pv.vec <- matrix(pv.de2, ncol = 1)
      qv.vec <- matrix(p.adjust(pv.vec[, 1], "fdr"), ncol = 1)
      fc.vec <- matrix(fc.de2, ncol = 1)
      m.vec <- matrix(rowMeans(prop0), ncol = 1)
      rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(m.vec) <- rownames(prop)
      colnames(pv.vec) <- "Pvalue"
      colnames(qv.vec) <- "Qvalue"
      colnames(fc.vec) <- "SpearmanCorr"
      colnames(m.vec) <- "Mean"
      pv.list[[LOI]] <- pv.vec
      qv.list[[LOI]] <- qv.vec
      fc.list[[LOI]] <- fc.vec
      m.list[[LOI]] <- m.vec
      if (method %in% c("perm", "perm.pair")) {
        qv.perm.list[[LOI]] <- t(t(perm.obj$p.adj.fdr))
        R2.list[[LOI]] <- t(t(perm.obj$R2))
        coef.list[[LOI]] <- perm.obj$coefficients
      }
      if (method %in% c("perm", "perm.pair")) {
        coef.mat <- perm.obj$coefficients
        colnames(coef.mat) <- paste0("coef_", colnames(coef.mat))
        res <- cbind(pv.vec, qv.vec, Qvalue.perm = perm.obj$p.adj.fdr,
                     R2 = perm.obj$R2, coef.mat, fc.vec)
        rownames(res) <- rownames(prop)
      } else {
        res <- cbind(m.vec, fc.vec, pv.vec, qv.vec)
        rownames(res) <- rownames(prop)
      }
      write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_",
                            LOI, "_", ann, ".csv"))
      if (mt.method == "fdr") {
        res.final <- rbind(res.final, res[res[, "Qvalue"] <= cutoff,
                                          , drop = F])
      }
      if (mt.method == "raw") {
        res.final <- rbind(res.final, res[res[, "Pvalue"] <= cutoff,
                                          , drop = F])
      }
      if (mt.method == "fdr.perm") {
        res.final <- rbind(res.final, res[res[, "Qvalue.perm"] <=
                                            cutoff, , drop = F])
      }
    }
    if (!is.null(res.final)) {
      colnames(res.final) <- colnames(res)
      write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_",
                                  mt.method, "_", cutoff, "_", ann, ".csv"))
    }
    return(list(pv.list = pv.list, fc.list = fc.list, qv.list = qv.list,
                m.list = m.list, qv.perm.list = qv.perm.list, coef.list = coef.list,
                R2.list = R2.list))
  }
}
confint.lme <- function(object, parm, level = 0.95, ...) {
  cf <- fixed.effects(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames else if (is.numeric(parm))
      parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- stats:::format.perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- sqrt(diag(vcov(object)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}
build.decision.tree <- function(data.obj, resp.name, taxa.level = "Species",
                                binary = FALSE, taxa, aug.var = NULL, ann = "All") {
  ann <- paste(taxa.level, ann, sep = "_")
  response <- data.obj$meta.dat[, resp.name]
  ct <- data.obj$abund.list[[taxa.level]]
  prop <- t(t(ct)/colSums(ct))
  prop <- prop[taxa, , drop = F]
  if (binary == TRUE) {
    prop <- (prop != 0)
  }
  dat <- as.data.frame(t(prop))
  # Rev: 2017_02_17 Add additional variables from meta dat
  if (!is.null(aug.var)) {
    dat <- cbind(dat, data.obj$meta.dat[, aug.var])
  }
  dat <- data.frame(dat, response)
  try(if (is.factor(response)) {
    fit <- rpart(response ~ ., method = "class", data = dat)
    post(fit, file = paste0("Taxa_Unpruned_Classification_tree_", ann,
                            ".ps"), title = "Unpruned Classification Tree")
    pfit <- prune(fit, cp = fit$cptable[which.min(fit$cptable[, "xerror"]),
                                        "CP"])
    post(pfit, file = paste0("Taxa_Pruned_Classification_", ann, ".ps"),
         title = "Pruned Classification Tree")
  } else {
    fit <- rpart(response ~ ., method = "anova", data = dat)
    post(fit, file = paste0("Taxa_Unpruned_Regression_tree_", ann,
                            ".ps"), title = "Unpruned Regression Tree")
    pfit <- prune(fit, cp = fit$cptable[which.min(fit$cptable[, "xerror"]),
                                        "CP"])
    post(pfit, file = paste0("Taxa_Pruned_Regression_tree_", ann, ".ps"),
         title = "Pruned RegressionTree")
  })
}
bootstrap.pwr <- function(pcs, formula, dat, ns = NULL, perm.no = 199,
                          iter.no = 100) {
  if (is.null(ns)) {
    ns <- nrow(dat) * seq(1, 4, len = 8)
  }
  pvs <- numeric(length(ns))
  names(pvs) <- paste(ns)
  for (n in ns) {
    cat(".")
    temp <- sapply(1:iter.no, function(i) {
      bt.ind <- sample(1:nrow(dat), n, repl = T)
      dat.bt <- dat[bt.ind, ]
      dist.bt <- dist(pcs[bt.ind, ])
      aov.tab <- adonis(as.formula(paste("dist.bt", formula)), dat = dat.bt,
                        permutations = perm.no)$aov.tab
      pv <- aov.tab[nrow(aov.tab) - 2, ncol(aov.tab)]
    })
    pvs[paste(n)] <- mean(temp <= 0.05)
  }
  return(pvs)
}
perform_power_analysis <- function(data.obj, dist.obj, dist.names = c("UniFrac",
                                                                      "GUniFrac", "WUniFrac", "BC"), formula = NULL, grp.name = NULL, adj.name = NULL,
                                   ann = "", ...) {
  if (is.null(formula)) {
    if (is.null(adj.name)) {
      formula <- paste("~", grp.name)
    } else {
      formula <- paste("~", paste(adj.name, collapse = "+"), "+",
                       grp.name)
    }
  }
  df <- data.obj$meta.dat
  pvm <- NULL
  for (dist.name in dist.names) {
    cat("*")
    dist.mat <- dist.obj[[dist.name]]
    pcs <- cmdscale(dist.mat, k = nrow(dist.mat) - 1)
    pvs <- bootstrap.pwr(pcs, formula, df, ...)
    pvm <- rbind(pvm, pvs)
  }
  colnames(pvm) <- names(pvs)
  rownames(pvm) <- dist.names
  pvdf <- melt(pvm)
  colnames(pvdf) <- c("Distance_type", "Sample_size", "Value")
  pdf(paste0("Power_curve_bootstrap_", ann, ".pdf"))
  g.obj <- ggplot(pvdf, aes(x = Sample_size, y = Value)) + geom_point() +
    geom_line() + ylab("Power") + facet_wrap(~Distance_type, ncol = 2) +
    theme_bw()
  print(g.obj)
  dev.off()
  return(pvdf)
}
# More robust identification
k.nonoverlap.se <- function(tab, SE.factor = 1) {
  max.v <- which.max(tab[, "gap"])
  for (i in (max.v - 1):1) {
    if (tab[i, "gap"] + SE.factor * tab[i, "SE.sim"] < tab[i + 1, "gap"] -
        SE.factor * tab[i + 1, "SE.sim"]) {
      break
    }
  }
  return(i + 1)
}
##########
B2M <- function(x) {
  x[x == 0] <- min(x[x != 0])
  x[x == 1] <- max(x[x != 1])
  log(x/(1 - x))
}
fastDist <- function(X) {
  temp <- colSums(X^2)
  D <- outer(temp, temp, "+") - 2 * t(X) %*% X
  diag(D) <- 0
  sqrt(D)
}
fastLM <- function(Y, M) {
  Y <- as.matrix(Y)
  XXI <- solve(t(M) %*% M)
  dof <- ncol(Y) - ncol(M)
  est <- XXI %*% t(M) %*% t(Y)
  resid <- t(Y) - M %*% est
  sigma <- sqrt(colSums(resid^2)/dof)
  Pvec <- 2 * pt(-abs(t(est/(sqrt(diag(XXI))))/sigma), dof)
  return(Pvec)
}
matrix.paste0 <- function(m.list) {
  p <- length(m.list)
  m <- max(unlist(sapply(m.list, function(x) nrow(x))))
  n <- max(unlist(sapply(m.list, function(x) ncol(x))))
  res <- NULL
  for (i in 1:p) {
    mat <- m.list[[i]]
    if (!is.matrix(mat)) {
      mat <- matrix(mat, m, n)
    } else {
      if (!is.null(rownames(mat))) {
        row.names <- rownames(mat)
      }
      if (!is.null(colnames(mat))) {
        col.names <- colnames(mat)
      }
    }
    res <- paste0(res, mat)
  }
  return(matrix(res, m, n, dimnames = list(row.names, col.names)))
}
# Typer I error 1
calculateK <- function(G) {
  n <- nrow(G)
  m <- mean(diag(G))
  v <- var(G[lower.tri(G)])
  sigma <- v/m
  k <- m/sigma
  se <- 1.96 * sqrt(2/n)
  list(K = k, K.lower = (sqrt(k) - se)^2, K.upper = (sqrt(k) + se)^2,
       sigma = sigma)
}
fPERMANOVA <- function(formula, data = NULL, contr.unordered = "contr.sum",
                       contr.ordered = "contr.poly") {
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  # The first group includes the intercept, nterms count the intercept
  u.grps <- unique(grps)
  nterms <- length(u.grps)
  Z <- rhs[, grps %in% u.grps[1:(nterms - 1)], drop = F]
  XZ <- rhs[, grps %in% u.grps[1:(nterms)], drop = F]
  n <- nrow(XZ)
  HZ <- Z %*% solve(t(Z) %*% (Z)) %*% t(Z)
  HXZ <- XZ %*% solve(t(XZ) %*% (XZ)) %*% t(XZ)
  HX <- HXZ - HZ
  HIX <- diag(n) - HXZ
  if (class(lhs) == "dist") {
    D <- -as.matrix(lhs)^2/2
  } else {
    D <- -lhs^2/2
  }
  G <- mean(D) + D - rowMeans(D) - matrix(rep(1, n), ncol = 1) %*% colMeans(D)
  G <- G/n
  obj <- calculateK(G)
  K <- obj$K
  sigma <- obj$sigma
  df1 <- ncol(XZ) - ncol(Z)
  df2 <- n - ncol(XZ)
  MSS <- sum(G * HX)
  RSS <- sum(G * HIX)
  TSS <- sum(diag(G))
  # RSS <- TSS - MSS
  f.stat <- (MSS/df1)/(RSS/df2)
  # p.value1 <- pnorm(f.stat, mean=1, sd=sqrt(2 / (K * df1)),
  # lower.tail=F)
  p.value2 <- pchisq(f.stat * K * df1, df = K * df1, lower.tail = F)
  SumsOfSqs <- c(MSS, RSS, TSS)
  tab <- data.frame(Df = c(df1, df2, n - 1), SumsOfSqs = n * SumsOfSqs,
                    MeanSqs = c(n * MSS/df1, n * RSS/df2, NA), F.Model = c(f.stat,
                                                                           NA, NA), R2 = c(MSS/TSS, NA, NA), `Pr(>F)` = c(p.value2, NA,
                                                                                                                          NA))
  rownames(tab) <- c("Model(Adjusted)", "Residuals", "Total")
  colnames(tab)[ncol(tab)] <- c("Pr(>F)")
  class(tab) <- c("anova", class(tab))
  attr(tab, "heading") <- c("F stat and P value of the last term is adjusted by previous terms!\n")
  out <- list(aov.tab = tab, K = K, sigma = sigma, call = match.call())
  out
}
# Determine the optimal depth to filter the samples
optim_dep <- function(data.obj, ID.ctrl, ID.other = NULL, taxa.level = "Species",
                      ann = "") {
  ID <- rownames(data.obj$meta.dat)
  if (is.null(ID.other)) {
    ID.other <- setdiff(ID, ID.ctrl)
  }
  dep <- colSums(data.obj$otu.tab)
  gen <- data.obj$abund.list[[taxa.level]]
  gen <- t(gen)/colSums(gen)
  dist.mat <- as.matrix(vegdist(gen))
  dist.mat <- dist.mat[ID.other, ID.ctrl]
  dep <- dep[ID.other]
  ord <- order(dep)
  dep <- dep[ord]
  dist.mat <- dist.mat[ord, ]
  # data <- data.frame(distance=as.vector(t(dist.mat)),
  # depth=rep(log10(dep), each = ncol(dist.mat)))
  data <- data.frame(distance = rowMeans(dist.mat), depth = log10(dep))
  pdf(paste0("Distance2CtrlvsDepth.", taxa.level, ".", ann, ".pdf"),
      width = 6, height = 5)
  obj <- ggplot(data, aes(x = depth, y = distance)) + # geom_boxplot() +
    geom_smooth() + geom_jitter() + geom_vline(xintercept = log10(2000),
                                               col = "red") + xlab("Log10(depth)") + theme_bw()
  print(obj)
  dev.off()
}
# Determine the influential taxa - Need to test
contaminant_detection <- function(data.obj, ID.ctrl, ID.other = NULL, taxa.level = "Species",
                                  ann = "") {
  ID <- rownames(data.obj$meta.dat)
  if (is.null(ID.other)) {
    ID.other <- setdiff(ID, ID.ctrl)
  }
  dep <- colSums(data.obj$otu.tab)
  gen <- data.obj$abund.list[[taxa.level]]
  gen <- t(gen)/colSums(gen)
  func <- function(gen, ID.other, ID.ctrl) {
    n.ctrl <- length(ID.ctrl)
    n.other <- length(ID.other)
    dist.mat <- as.matrix(vegdist(gen))
    dist1 <- rowMeans(dist.mat[ID.other, ID.ctrl])
    # dist2 <- rowMeans(dist.mat[ID.other, ID.other]) * n.other / (n.other
    # - 1)
    # mean(dist1) / mean(dist2)
    mean(dist1)
  }
  r0 <- func(gen, ID.other, ID.ctrl)
  r.vec <- numeric(ncol(gen))
  names(r.vec) <- colnames(gen)
  for (i in 1:ncol(gen)) {
    if (i%%10 == 0)
      cat(".")
    # gen2 <- gen[, -i] gen2 <- gen2 / rowSums(gen2)
    gen2 <- gen
    gen2[, i] <- sample(gen2[, i])
    r.vec[i] <- func(gen2, ID.other, ID.ctrl)
  }
  return(list(r0 = r0, r.vec = r.vec))
}
