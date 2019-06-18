perform_cluster_analysis <- function(data.obj, dist.obj, dist.name = c("UniFrac"),
                                     k.best = NULL, method = "pam", stat = "gap", grp.name = NULL, adj.name = NULL,
                                     subject = NULL, ann = "", seed = 1234) {
  df <- data.obj$meta.dat
  if (!is.null(grp.name)) {
    grp <- df[, grp.name]
  }
  if (!is.null(adj.name)) {
    adj <- df[, adj.name]
  }
  cat(dist.name, " distance ...\n")
  mat <- dist.obj[[dist.name]]
  if (is.null(k.best)) {
    pc.obj <- cmdscale(as.dist(mat), k = ncol(mat) - 1, eig = T)
    eig <- pc.obj$eig
    eig <- eig[eig > 0]
    pvar <- eig/sum(eig)
    pc <- pc.obj$points
    cat("Assess cluster number by gap statistics ...\n")
    set.seed(seed)
    gs <- gapstat_ord(pc, axes = 1:ncol(pc), verbose = FALSE)
    plot_clusgap(gs)
    ggsave(paste0("cluster_assess_gap_statistic_", dist.name, ann,
                  ".pdf"), width = 6, height = 6)
    # print(gs, method='Tibs2001SEmax') print(gs, method='firstSEmax')
    tab <- gs$Tab
    write.csv(tab, paste0("gap_stat_", dist.name, ".csv"))
    k.best.gap <- k.nonoverlap.se(tab)
    cat("Gap statistic finds ", k.best.gap, " clusters.\n")
    cat("Assess cluster number by average silhouette width ...\n")
    asw <- numeric(20)
    for (k in 2:20) {
      asw[k] <- pam(as.dist(mat), k)$silinfo$avg.width
    }
    k.best.asw <- which.max(asw)
    cat("silhouette-optimal number of clusters:", k.best.asw, "\n")
    pdf(paste0("cluster_assess_asw_statistic_", dist.name, ann, ".pdf"),
        width = 6, height = 6)
    plot(1:20, asw, type = "h", main = "pam() clustering assessment",
         xlab = "k (# clusters)", ylab = "average silhouette width")
    axis(1, k.best.asw, paste("best", k.best.asw, sep = "\n"), col = "red",
         col.axis = "red")
    dev.off()
    cat("ASW statistic finds ", k.best.asw, " clusters.\n")
    if (stat == "gap") {
      k.best <- k.best.gap
    } else {
      if (stat == "asw") {
        k.best <- k.best.asw
      }
    }
  }
  if (k.best == 1) {
    cat("No robust cluster structure found!\n")
  } else {
    pam.obj <- pam(as.dist(mat), k = k.best)
    pam.class <- factor(pam.obj$clustering)
    write.csv(pam.class, paste0("Cluster.membership", dist.name, ".csv"))
    df$Cluster <- pam.class
    data.obj$meta.dat <- df
    if (is.null(grp.name)) {
      generate_ordination(data.obj, dist.obj, dist.name, grp.name = "Cluster",
                          ann = paste0("Cluster", dist.name, ann))
    } else {
      generate_ordination(data.obj, dist.obj, dist.name, grp.name = "Cluster",
                          strata = grp.name, ann = paste0("Cluster", dist.name, ann))
    }
    # Define cluster characteristics
    diff.obj <- perform_differential_analysis(data.obj, grp.name = "Cluster",
                                              taxa.levels = c("Genus"), method = "kruskal", mt.method = "fdr",
                                              cutoff = 0.01, prev = 0.1, minp = 0.002, ann = "Cluster")
    visualize_differential_analysis(data.obj, diff.obj, grp.name = "Cluster",
                                    taxa.levels = c("Genus"), indivplot = FALSE, mt.method = "fdr",
                                    cutoff = 0.01, ann = "Cluster")
    try(if (!is.null(grp.name)) {
      if (is.null(adj.name)) {
        form <- as.formula(paste("yy ~", grp.name))
      } else {
        form <- as.formula(paste("yy ~", paste(adj.name, collapse = "+"),
                                 "+", grp.name))
      }
      # Enrichment analysis - logistic regression - 1 vs other
      cat("Testing for association with the clusters - 1 vs other ...\n")
      sink(paste0("Cluster_association_test", dist.name, ann, ".txt"))
      if (is.null(subject)) {
        cat("Generalized linear model (logistic regression) is performed.\n")
        for (clus in levels(df$Cluster)[1:(nlevels(df$Cluster))]) {
          cat("Test for enrichment in cluster", clus, "\n")
          y <- rep(0, length(df$Cluster))
          y[df$Cluster == clus] <- 1
          df$yy <- y
          prmatrix(summary(glm(form, data = df, family = binomial))$coefficients)
        }
      } else {
        cat("Generalized linear mixed effects model (logistic regression) is performed.\n")
        for (clus in levels(df$Cluster)[1:(nlevels(df$Cluster))]) {
          cat("Test for enrichment in cluster", clus, "\n")
          y <- rep(0, length(df$Cluster))
          y[df$Cluster == clus] <- 1
          df$yy <- y
          prmatrix(summary(glmmPQL(form, data = df, random = as.formula(paste0("~ 1|",
                                                                               subject)), family = binomial, verbose = F))$tTable)
        }
      }
      sink()
    })
  }
}
