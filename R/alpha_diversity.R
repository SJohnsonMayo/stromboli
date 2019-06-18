# Rev: 2016_12_25 Add anova, and record the results Rev: 2018_02_24
# Remove dependence on Phyloseq
perform_alpha_test <- function(data.obj, phylo.obj = NULL, rarefy = TRUE,
                               depth = NULL, iter.no = 5, measures = c("Observed", "Chao1", "Shannon",
                                                                       "InvSimpson"), model = "lm", formula = NULL, grp.name = NULL, adj.name = NULL,
                               ann = "", seed = 123, ...) {
  # Implement future random effects model Rev: 2017_08_23
  if (is.null(phylo.obj)) {
    phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows = T),
                          phy_tree(data.obj$tree), tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
  }
  result <- list()
  if (rarefy == TRUE) {
    if (is.null(depth)) {
      depth <- min(sample_sums(phylo.obj))
    } else {
      if (depth > min(sample_sums(phylo.obj))) {
        ind <- sample_sums(phylo.obj) >= depth
        cat(sum(!ind), " samples do not have sufficient number of reads!\n")
        sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ]
        data.obj <- subset_data(data.obj, ind)
      }
    }
    x <- 0
    sink("temp.txt")
    for (i in 1:iter.no) {
      phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed = 12345 +
                                        i)
      x <- x + estimate_richness(phylo.even, measures = measures)
    }
    sink()
    x <- x/iter.no
  } else {
    x <- estimate_richness(phylo.obj, measures = measures)
  }
  df <- data.obj$meta.dat
  fitted.obj <- list()
  if (rarefy == T) {
    sink(paste0("Alpha_diversity_test_results_rarefied_", ann, ".txt"))
  } else {
    sink(paste0("Alpha_diversity_test_results_unrarefied_", ann, ".txt"))
  }
  date()
  # variables to adjust always come first in anova analyses
  if (is.null(formula)) {
    if (is.null(adj.name)) {
      formula <- paste("~", grp.name)
    } else {
      formula <- paste("~", paste(adj.name, collapse = "+"), "+",
                       grp.name)
    }
  }
  for (measure in measures) {
    cat("Alpha diversity:", measure, "\n")
    xx <- x[, measure]
    if (model == "lm") {
      cat("Linear model:\n")
      lm.obj <- lm(as.formula(paste("xx ", formula)), df, ...)
      prmatrix(summary(lm.obj)$coefficients)
      cat("\nANOVA:\n")
      print(anova(lm.obj))
      cat("\n")
      fitted.obj[[measure]] <- lm.obj
    }
    if (model == "lme") {
      df$xx <- xx
      cat("Linear mixed effects model:\n")
      lm.obj <- lme(as.formula(paste("xx ", formula)), df, method = "ML",
                    ...)
      prmatrix(summary(lm.obj)$tTable)
      cat("\nANOVA:\n")
      print(anova(lm.obj))
      cat("\n")
      fitted.obj[[measure]] <- lm.obj
    }
    cat("\n")
  }
  sink()
  result$fitted.obj <- fitted.obj
  result$alpha.diversity <- x
  # Rev: 2016_12_25
  return(invisible(result))
}
# New: 2018_02_25
generate_alpha_diversity <- function(data.obj, rarefy = TRUE, depth = NULL,
                                     iter.no = 5, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
                                     seed = 123) {
  # Rev: 2017_08_23
  OTU <- t(data.obj$otu.tab)
  depths <- rowSums(OTU)
  result <- list()
  if (rarefy == TRUE) {
    if (is.null(depth)) {
      depth <- min(depths)
    } else {
      if (depth > min(depths)) {
        ind <- depths >= depth
        cat(sum(!ind), " samples do not have sufficient number of reads!\n")
        # data.obj <- subset_data(data.obj, ind)
      }
    }
    set.seed(123)
    x <- 0
    for (i in 1:iter.no) {
      OTU2 <- Rarefy(OTU, depth = depth)$otu.tab.rff
      x <- x + estimate_richness_(OTU2, measures = measures)
    }
    x <- x/iter.no
    rownames(x) <- rownames(OTU2)
  } else {
    x <- estimate_richness_(OTU, measures = measures)
    rownames(x) <- rownames(OTU)
  }
  return(x)
}
# Rev: 2016_09_10 Rev: 2016_11_28 Rev: 2017_04_18 Rev: 2018_03_06
# position_jitter(h=0)
generate_alpha_boxplot <- function(data.obj, phylo.obj = NULL, rarefy = TRUE,
                                   depth = NULL, grp.name, strata = NULL, measures = c("Observed", "Chao1",
                                                                                       "Shannon", "InvSimpson"), gg.cmd = NULL, ann = "", subject = NULL,
                                   p.size = 2.5, l.size = 0.5, hei = NULL, wid = NULL) {
  # Rev: 2017_08_23
  if (is.null(phylo.obj)) {
    phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows = T),
                          phy_tree(data.obj$tree), tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
  }
  # To be completed - jetter when strata is not null
  if (rarefy == TRUE) {
    if (is.null(depth)) {
      depth <- min(sample_sums(phylo.obj))
    } else {
      if (depth > min(sample_sums(phylo.obj))) {
        ind <- (sample_sums(phylo.obj) >= depth)
        cat(sum(!ind), " samples do not have sufficient number of reads!\n")
        sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ]
        data.obj <- subset_data(data.obj, ind)
      }
    }
    phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed = 12345)
    x <- estimate_richness(phylo.even, measures = measures)
  } else {
    x <- estimate_richness(phylo.obj, measures = measures)
  }
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  if (!is.null(subject)) {
    ID <- df[, subject]
  }
  if (is.null(hei) | is.null(wid)) {
    hei <- 5
    if (is.null(strata)) {
      wid <- 5
    } else {
      wid <- 5.5
    }
  }
  if (rarefy == T) {
    pdf(paste0("Alpha_diversity_boxplot_rarefied_", ann, ".pdf"), height = hei,
        width = wid)
  } else {
    pdf(paste0("Alpha_diversity_boxplot_unrarefied_", ann, ".pdf"),
        height = hei, width = wid)
  }
  if (is.null(subject)) {
    if (is.null(strata)) {
      for (measure in measures) {
        cat(measure, "\n")
        xx <- x[, measure]
        df2 <- data.frame(Value = xx, Group = grp)
        dodge <- position_dodge(width = 0.75)
        obj <- ggplot(df2, aes(x = Group, y = Value, col = Group)) +
          geom_boxplot(position = dodge, outlier.colour = NA) +
          geom_jitter(alpha = 0.6, size = 3, position = position_jitter(w = 0.1,
                                                                        h = 0)) + labs(y = measure) + theme(legend.position = "none")
        if (!is.null(gg.cmd)) {
          obj <- obj + eval(parse(text = gg.cmd))
        }
        print(obj)
      }
    } else {
      for (measure in measures) {
        cat(measure, "\n")
        xx <- x[, measure]
        grp2 <- df[, strata]
        df2 <- data.frame(Value = xx, Group = grp, Strata = grp2)
        dodge <- position_dodge(width = 0.9)
        obj <- ggplot(df2, aes(x = Strata, y = Value, col = Group,
                               fill = Group)) + geom_jitter(alpha = 0.6, size = 3, position = position_jitterdodge(dodge.width = 0.9)) +
          geom_boxplot(position = dodge, outlier.colour = NA, fill = "white",
                       alpha = 0.6) + labs(y = measure, x = strata)
        if (!is.null(gg.cmd)) {
          obj <- obj + eval(parse(text = gg.cmd))
        }
        print(obj)
      }
    }
    dev.off()
  } else {
    if (is.null(strata)) {
      for (measure in measures) {
        cat(measure, "\n")
        xx <- x[, measure]
        df2 <- data.frame(Value = xx, Group = grp, subject = ID)
        dodge <- position_dodge(width = 0.75)
        obj <- ggplot(df2, aes(x = Group, y = Value, shape = Group,
                               group = subject)) + geom_point(size = p.size) + geom_line(size = l.size) +
          labs(y = measure) + theme(legend.position = "none")
        if (!is.null(gg.cmd)) {
          obj <- obj + eval(parse(text = gg.cmd))
        }
        print(obj)
      }
    } else {
      warning("Stratum is disallowed when subject is specified!\n")
    }
    dev.off()
  }
}
# Rev: 2018_03_09 Removing phyloseq dependency and add 'subject'
# parameter, and remove 'model' parameter
generate_rarefy_curve2 <- function(data.obj, grp.name, measures = c("Observed",
                                                                    "Chao1", "Shannon", "InvSimpson"), depth = NULL, iter.no = 5, npoint = 10,
                                   seed = 123, ann = "", gg.cmd = "theme(legend.justification=c(1,0), legend.position=c(1,0))",
                                   wid = 5, hei = 5, error = "se") {
  cat("Create rarefaction curves!\n")
  if (is.null(depth)) {
    depth <- min(colSums(data.obj$otu.tab))
  }
  if (error == "se")
    k <- 1
  if (error == "ci")
    k <- 1.96
  ind <- colSums(data.obj$otu.tab) >= depth
  if (sum(!ind) != 0) {
    cat(sum(!ind), " samples do not have sufficient number of reads!\n")
    data.obj <- subset_data(data.obj, ind)
  }
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  if (is.character(grp)) {
    grp <- factor(grp)
  }
  if (!is.factor(grp)) {
    stop("Rarefaction curve needs a factor!\n")
  }
  res <- NULL
  incr <- depth%/%npoint
  for (dep in c(10, incr * (1:npoint))) {
    x <- generate_alpha_diversity(data.obj, rarefy = TRUE, depth = dep,
                                  iter.no = iter.no, measures = measures, seed = seed)
    res <- rbind(res, t(x[, measures, drop = F]))
  }
  colnames(res) <- rownames(df)
  cat("Wait ...\n")
  pdf(paste0("Alpha_diversity_Rarefaction_Curve_", ann, ".pdf"), width = wid,
      height = hei)
  for (i in 1:length(measures)) {
    measure <- measures[i]
    cat("Measure: ", measure, "\n")
    res2 <- res[(0:(npoint)) * length(measures) + i, , drop = F]
    m <- t(apply(res2, 1, function(x) tapply(x, grp, mean)))
    se <- t(apply(res2, 1, function(x) tapply(x, grp, function(y) k *
                                                sd(y)/sqrt(length(y)))))
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
                    position = position_dodge(0.2)) + geom_line() + geom_point(size = 3,
                                                                               shape = 21, fill = "white") + labs(y = measure)
    if (!is.null(gg.cmd)) {
      obj <- obj + eval(parse(text = gg.cmd))
    }
    print(obj)
  }
  dev.off()
  cat("Finished!\n")
}
# New: 2018_03_09 Removing phyloseq dependency and add 'subject'
# parameter, and remove 'model' parameter
perform_alpha_test2 <- function(data.obj, alpha.obj = NULL, rarefy = TRUE,
                                depth = NULL, iter.no = 5, measures = c("Observed", "Chao1", "Shannon",
                                                                        "InvSimpson"), model = "lm", formula = NULL, grp.name = NULL, adj.name = NULL,
                                subject = NULL, ann = "", seed = 123, ...) {
  if (is.null(alpha.obj)) {
    alpha.obj <- generate_alpha_diversity(data.obj, rarefy = rarefy,
                                          depth = depth, iter.no = iter.no, measures = measures, seed = seed)
  } else {
    if (sum(!(rownames(alpha.obj) %in% rownames(data.obj$meta.dat))) !=
        0) {
      stop("alpha.obj contains samples not in data.obj!\n")
    }
  }
  x <- alpha.obj
  df <- data.obj$meta.dat[rownames(alpha.obj), ]
  if (is.null(subject)) {
    model <- "lm"
  } else {
    model <- "lme"
  }
  result <- list()
  fitted.obj <- list()
  if (rarefy == T) {
    sink(paste0("Alpha_diversity_test_results_rarefied_", ann, ".txt"))
  } else {
    sink(paste0("Alpha_diversity_test_results_unrarefied_", ann, ".txt"))
  }
  date()
  # variables to adjust always come first in anova analyses
  if (is.null(formula)) {
    if (is.null(adj.name)) {
      formula <- paste("~", grp.name)
    } else {
      formula <- paste("~", paste(adj.name, collapse = "+"), "+",
                       grp.name)
    }
  }
  for (measure in measures) {
    cat("Alpha diversity:", measure, "\n")
    xx <- x[, measure]
    if (model == "lm") {
      cat("Linear model:\n")
      lm.obj <- lm(as.formula(paste("xx ", formula)), df, ...)
      prmatrix(summary(lm.obj)$coefficients)
      cat("\nANOVA:\n")
      print(anova(lm.obj))
      cat("\n")
      fitted.obj[[measure]] <- lm.obj
    }
    if (model == "lme") {
      df$xx <- xx
      cat("Linear mixed effects model:\n")
      lm.obj <- lme(as.formula(paste("xx ", formula)), df, method = "ML",
                    random = as.formula(paste0(" ~ 1 | ", subject)), ...)
      prmatrix(summary(lm.obj)$tTable)
      cat("\nANOVA:\n")
      print(anova(lm.obj))
      cat("\n")
      fitted.obj[[measure]] <- lm.obj
    }
    cat("\n")
  }
  sink()
  result$fitted.obj <- fitted.obj
  result$alpha.diversity <- x
  # Rev: 2016_12_25
  return(invisible(result))
}
# New: 2018_03_09 Removing phyloseq dependency
generate_alpha_boxplot2 <- function(data.obj, alpha.obj = NULL, rarefy = TRUE,
                                    depth = NULL, iter.no = 5, measures = c("Observed", "Chao1", "Shannon",
                                                                            "InvSimpson"), seed = 123, grp.name, strata = NULL, gg.cmd = NULL,
                                    ann = "", subject = NULL, p.size = 2.5, l.size = 0.5, hei = NULL, wid = NULL) {
  # Rev: 2017_08_23
  if (is.null(alpha.obj)) {
    alpha.obj <- generate_alpha_diversity(data.obj, rarefy = rarefy,
                                          depth = depth, iter.no = iter.no, measures = measures, seed = seed)
  } else {
    if (sum(!(rownames(alpha.obj) %in% rownames(data.obj$meta.dat))) !=
        0) {
      stop("alpha.obj contains samples not in data.obj!\n")
    }
  }
  x <- alpha.obj
  df <- data.obj$meta.dat[rownames(alpha.obj), ]
  grp <- df[, grp.name]
  if (!is.null(subject)) {
    ID <- df[, subject]
  }
  if (is.null(hei) | is.null(wid)) {
    hei <- 5
    if (is.null(strata)) {
      wid <- 5
    } else {
      wid <- 5.5
    }
  }
  if (rarefy == T) {
    pdf(paste0("Alpha_diversity_boxplot_rarefied_", ann, ".pdf"), height = hei,
        width = wid)
  } else {
    pdf(paste0("Alpha_diversity_boxplot_unrarefied_", ann, ".pdf"),
        height = hei, width = wid)
  }
  if (is.null(subject)) {
    if (is.null(strata)) {
      for (measure in measures) {
        cat(measure, "\n")
        xx <- x[, measure]
        df2 <- data.frame(Value = xx, Group = grp)
        dodge <- position_dodge(width = 0.75)
        obj <- ggplot(df2, aes(x = Group, y = Value, col = Group)) +
          geom_boxplot(position = dodge, outlier.colour = NA, lwd = 0.35) +
          geom_jitter(alpha = 0.6, size = p.size, position = position_jitter(w = 0.1,
                                                                             h = 0)) + labs(y = measure) + theme(legend.position = "none")
        if (!is.null(gg.cmd)) {
          obj <- obj + eval(parse(text = gg.cmd))
        }
        print(obj)
      }
    } else {
      for (measure in measures) {
        cat(measure, "\n")
        xx <- x[, measure]
        grp2 <- df[, strata]
        df2 <- data.frame(Value = xx, Group = grp, Strata = grp2)
        dodge <- position_dodge(width = 0.75)
        obj <- ggplot(df2, aes(x = Strata, y = Value, col = Group,
                               fill = Group)) + geom_boxplot(position = dodge, outlier.colour = NA,
                                                             fill = "white", alpha = 0.65, lwd = 0.35, width = 0.65) +
          geom_point(position = position_jitterdodge(dodge.width = 0.75),
                     size = p.size, alpha = 0.6) + labs(y = measure, x = strata)
        if (!is.null(gg.cmd)) {
          obj <- obj + eval(parse(text = gg.cmd))
        }
        print(obj)
      }
    }
    dev.off()
  } else {
    if (is.null(strata)) {
      for (measure in measures) {
        cat(measure, "\n")
        xx <- x[, measure]
        df2 <- data.frame(Value = xx, Group = grp, subject = ID)
        dodge <- position_dodge(width = 0.75)
        obj <- ggplot(df2, aes(x = Group, y = Value, shape = Group,
                               group = subject)) + geom_point(size = p.size) + geom_line(size = l.size) +
          labs(y = measure) + theme(legend.position = "none")
        if (!is.null(gg.cmd)) {
          obj <- obj + eval(parse(text = gg.cmd))
        }
        print(obj)
      }
    } else {
      warning("Stratum is disallowed when subject is specified!\n")
    }
    dev.off()
  }
}
# New: 2018_03_08 Removing phyloseq dependency Rev: 2018_07_16 More
# advanced strata support
generate_alpha_scatterplot2 <- function(data.obj, alpha.obj = NULL, grp.name,
                                        subject = NULL, rarefy = TRUE, depth = NULL, iter.no = 5, seed = 123,
                                        measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), strata = NULL,
                                        strata2 = NULL, combine.strata = FALSE, combine.strata2 = FALSE, smooth.method = "loess",
                                        pt.shape = 16, pt.alpha = 0.5, pt.size = 2, subject.pt = FALSE, gg.cmd = NULL,
                                        ann = "", hei = NULL, wid = NULL, pdf = TRUE) {
  if (is.null(alpha.obj)) {
    alpha.obj <- generate_alpha_diversity(data.obj, rarefy = rarefy,
                                          depth = depth, iter.no = iter.no, measures = measures, seed = seed)
  } else {
    if (sum(!(rownames(alpha.obj) %in% rownames(data.obj$meta.dat))) !=
        0) {
      stop("alpha.obj contains samples not in data.obj!\n")
    }
  }
  x <- alpha.obj
  df <- data.obj$meta.dat[rownames(alpha.obj), ]
  grp <- df[, grp.name]
  if (is.null(subject)) {
    ID <- NA
  } else {
    ID <- factor(df[, subject])
  }
  if (is.null(hei) | is.null(wid)) {
    hei <- 5
    wid <- 6
    if (!is.null(strata) & combine.strata == FALSE) {
      wid <- 4 * nlevels(factor(df[, strata]))
    }
    if (!is.null(strata2) & combine.strata2 == FALSE) {
      hei <- 2.5 * nlevels(factor(df[, strata2]))
    }
  }
  if (pdf == TRUE) {
    if (rarefy == T) {
      pdf(paste0("Alpha_diversity_scatterplot_rarefied_", ann, ".pdf"),
          height = hei, width = wid)
    } else {
      pdf(paste0("Alpha_diversity_scatterplot_unrarefied_", ann,
                 ".pdf"), height = hei, width = wid)
    }
  }
  if (is.null(strata)) {
    for (measure in measures) {
      cat(measure, "\n")
      xx <- x[, measure]
      df2 <- data.frame(Value = xx, Group = grp, ID = ID)
      obj <- ggplot(df2, aes(x = Group, y = Value, group = ID)) +
        # geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
        geom_smooth(method = smooth.method, aes(group = NA), size = 0.75) +
        labs(y = "Value", x = grp.name, title = measure) + theme(legend.position = "none")
      if (!is.null(subject)) {
        if (subject.pt == FALSE) {
          obj <- obj + geom_line(size = 0.25, alpha = 0.75)
        } else {
          obj <- obj + geom_line(size = 0.25, alpha = 0.75) + geom_point(size = pt.size,
                                                                         shape = pt.shape, alpha = pt.alpha)
        }
      } else {
        obj <- obj + geom_point(size = pt.size, shape = pt.shape,
                                alpha = pt.alpha)
      }
      if (!is.null(gg.cmd)) {
        obj <- obj + eval(parse(text = gg.cmd))
      }
      print(obj)
    }
  } else {
    for (measure in measures) {
      cat(measure, "\n")
      xx <- x[, measure]
      grp2 <- factor(df[, strata])
      if (!is.null(strata2)) {
        grp3 <- factor(df[, strata2])
      } else {
        grp3 <- grp2
      }
      grp4 <- factor(paste(grp2, grp3))
      # VERY bad names with only difference in capital letters, variable
      # names are also confusing Rev: 2018_07_16
      df2 <- data.frame(Value = xx, Group = grp, Strata = grp2, Strata2 = grp3,
                        Strata3 = grp4, ID = ID)
      if (is.null(strata2)) {
        if (combine.strata == TRUE) {
          obj <- ggplot(df2, aes(x = Group, y = Value, group = ID,
                                 col = Strata)) + # geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
            geom_smooth(method = smooth.method, aes(group = Strata),
                        size = 0.75) + labs(y = "Value", x = grp.name, title = measure)
        } else {
          obj <- ggplot(df2, aes(x = Group, y = Value, group = ID)) +
            # geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
            geom_smooth(method = smooth.method, aes(group = NA),
                        size = 0.75) + labs(y = "Value", x = grp.name, title = measure) +
            facet_grid(. ~ Strata)
        }
      } else {
        if (combine.strata == TRUE) {
          if (combine.strata2 == TRUE) {
            obj <- ggplot(df2, aes(x = Group, y = Value, group = ID,
                                   col = Strata3)) + # geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
              geom_smooth(method = smooth.method, aes(group = Strata3),
                          size = 0.75) + labs(y = "Value", x = grp.name, title = measure)
          } else {
            obj <- ggplot(df2, aes(x = Group, y = Value, group = ID,
                                   col = Strata)) + # geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
              geom_smooth(method = smooth.method, aes(group = Strata),
                          size = 0.75) + labs(y = "Value", x = grp.name, title = measure) +
              facet_grid(Strata2 ~ .)
          }
        } else {
          if (combine.strata2 == TRUE) {
            warning("Currently, combine.strata2 will be forced to be FALSE if combine.strata=FALSE!\n")
          }
          obj <- ggplot(df2, aes(x = Group, y = Value, group = ID)) +
            # geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
            geom_smooth(method = smooth.method, aes(group = NA),
                        size = 0.75) + labs(y = "Value", x = grp.name, title = measure) +
            facet_grid(Strata2 ~ Strata)
        }
      }
      obj <- obj + theme(legend.title = element_blank())
      if (!is.null(subject)) {
        if (subject.pt == FALSE) {
          obj <- obj + geom_line(size = 0.25, alpha = 0.75)
        } else {
          obj <- obj + geom_line(size = 0.25, alpha = 0.75) + geom_point(size = pt.size,
                                                                         shape = pt.shape, alpha = pt.alpha)
        }
      } else {
        obj <- obj + geom_point(size = pt.size, shape = pt.shape,
                                alpha = pt.alpha)
      }
      if (!is.null(gg.cmd)) {
        obj <- obj + eval(parse(text = gg.cmd))
      }
      print(obj)
    }
  }
  if (pdf == TRUE) {
    dev.off()
  }
}
# New: 2017_02_22 Depend on phyloseq - outdated. Please see the newer
# version generate_alpha_scatterplot2
generate_alpha_scatterplot <- function(data.obj, phylo.obj = NULL, rarefy = TRUE,
                                       depth = NULL, grp.name, strata = NULL, strata2 = NULL, wrap.nrow = 2,
                                       smooth.method = "auto", hei0 = NULL, wid0 = NULL, measures = c("Observed",
                                                                                                      "Chao1", "Shannon", "InvSimpson"), gg.cmd = NULL, ann = "") {
  # Rev: 2017_08_23
  if (is.null(phylo.obj)) {
    phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows = T),
                          phy_tree(data.obj$tree), tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
  }
  if (rarefy == TRUE) {
    if (is.null(depth)) {
      depth <- min(sample_sums(phylo.obj))
    } else {
      if (depth > min(sample_sums(phylo.obj))) {
        ind <- (sample_sums(phylo.obj) >= depth)
        cat(sum(!ind), " samples do not have sufficient number of reads!\n")
        sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ]
        data.obj <- subset_data(data.obj, ind)
      }
    }
    phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed = 12345)
    x <- estimate_richness(phylo.even, measures = measures)
  } else {
    x <- estimate_richness(phylo.obj, measures = measures)
  }
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  if (is.null(hei0)) {
    hei <- 5
  } else {
    hei <- hei0
  }
  if (is.null(strata)) {
    if (is.null(wid0)) {
      wid <- 5
    } else {
      wid <- wid0
    }
  } else {
    if (is.null(wid0)) {
      wid <- 5.5
    } else {
      wid <- wid0
    }
  }
  if (rarefy == T) {
    pdf(paste0("Alpha_diversity_scatterplot_rarefied_", ann, ".pdf"),
        height = hei, width = wid)
  } else {
    pdf(paste0("Alpha_diversity_scatterplot_unrarefied_", ann, ".pdf"),
        height = hei, width = wid)
  }
  if (is.null(strata)) {
    for (measure in measures) {
      cat(measure, "\n")
      xx <- x[, measure]
      df2 <- data.frame(Value = xx, Group = grp)
      obj <- ggplot(df2, aes(x = Group, y = Value)) + geom_point() +
        geom_smooth(method = smooth.method) + labs(y = "Value",
                                                   x = grp.name, title = measure) + theme(legend.position = "none")
      if (!is.null(gg.cmd)) {
        obj <- obj + eval(parse(text = gg.cmd))
      }
      print(obj)
    }
  } else {
    for (measure in measures) {
      cat(measure, "\n")
      xx <- x[, measure]
      grp2 <- df[, strata]
      if (!is.null(strata2)) {
        grp3 <- df[, strata2]
      } else {
        grp3 <- grp2
      }
      df2 <- data.frame(Value = xx, Group = grp, Strata = grp2, Strata2 = grp3)
      if (is.null(strata2)) {
        obj <- ggplot(df2, aes(x = Group, y = Value)) + geom_point() +
          geom_smooth(method = smooth.method) + labs(y = "Value",
                                                     x = grp.name, title = measure) + facet_wrap(~Strata,
                                                                                                 nrow = wrap.nrow)
      } else {
        obj <- ggplot(df2, aes(x = Group, y = Value)) + geom_point() +
          geom_smooth(method = smooth.method) + labs(y = "Value",
                                                     x = grp.name, title = measure) + facet_grid(Strata2 ~
                                                                                                   Strata)
      }
      obj <- obj + theme(legend.title = element_blank())
      if (!is.null(gg.cmd)) {
        obj <- obj + eval(parse(text = gg.cmd))
      }
      print(obj)
    }
  }
  dev.off()
}
