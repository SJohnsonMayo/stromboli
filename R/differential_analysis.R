
perform_taxa_compare_test <- function (data.obj, grp.name, level1, level2, level3, alternative='greater', nperm=1000, seed=123,
                                       taxa.levels=c('Phylum', 'Family', 'Genus'), taxa.name='All', prev=0.1, minp=0.002, ann='') {

  df <- data.obj$meta.dat
  grp <- df[, grp.name]

  grp.levels <- levels(grp)
  grp.nlevels <- nlevels(grp)
  ind1 <- which(grp == level1)
  ind2 <- which(grp == level2)
  ind3 <- which(grp == level3)

  pv.list <- list()

  for (LOI in taxa.levels) {
    cat(LOI, "\n")
    prop <- data.obj$abund.list[[LOI]]
    prop <- t(t(prop) / colSums(prop))
    pv.vec <- NULL
    if (taxa.name == 'All') {
      prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
    } else {
      prop <- prop[taxa.name, , drop=FALSE]
    }

    for (taxon in rownames(prop)) {
      cat(".")
      dist.mat <- as.matrix(dist(rank(prop[taxon, ])))
      pv <- distance_compare_test(dist.mat, ind1, ind2, ind3, alternative, nperm)
      pv.vec <- c(pv.vec, pv)
    }
    names(pv.vec) <- rownames(prop)
    pv.list[[LOI]] <- pv.vec
    cat("\n")
  }

  for (LOI in taxa.levels) {
    pv.vec <- pv.list[[LOI]]

    qv.vec <- p.adjust(pv.vec, 'fdr')


    write.csv(data.frame("p value"=pv.vec, "q value"=qv.vec), paste0("Taxa_TrendAnalysis_", LOI, "_", ann, ".csv"))
  }

}

# New: 2017_11_03
# Now default is averaged distance for within or between-subject comparison
perform_site_correlation_test <- function (data.obj, dist.obj, site.name, sites, subject,
                                           dist.names = c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), dist.ave = TRUE, nPerm = 999, ann = '') {

  pmat <- array(NA, c(length(dist.names), length(sites), length(sites)), dimnames=list(dist.names, sites, sites))
  for (site1 in sites) {
    for(site2 in sites) {
      if (site1 != site2 & is.na(pmat[1, site1, site2])) {
        cat('.')
        sampling.site <- c(site1, site2)
        samIDs <- data.obj$meta.dat[, site.name] %in% sampling.site
        data.obj1 <- subset_data(data.obj, samIDs)
        dist.obj1 <- subset_dist(dist.obj, samIDs)

        # Compute the observed distance within subject
        lIDs <- as.character(data.obj1$meta.dat[, site.name])
        pIDs <- as.character(data.obj1$meta.dat[, subject])
        uIDs <- as.character(unique(data.obj1$meta.dat[, subject]))


        pdf(paste('Permutation_test', site1, site2, 'Correlation_', ann, '.pdf', sep='_'), width=6, height=6)
        for (dist.name in dist.names) {
          cat(dist.name, '...\n')

          dist.mat <- dist.obj1[[dist.name]]
          # To self
          dist.o1 <- NULL
          for (uID in uIDs) {
            ind1 <- pIDs == uID & lIDs == site1
            ind2 <- pIDs == uID & lIDs == site2

            if (sum(ind1) != 0 & sum(ind2) != 0)
              if (dist.ave) {
                dist.o1 <- c(dist.o1, mean(dist.mat[ind1, ind2]))
              } else {
                dist.o1 <- c(dist.o1, dist.mat[ind1, ind2])
              }

          }

          # To other people
          dist.o2 <- NULL
          for (uID in uIDs) {
            ind1 <- pIDs == uID & lIDs == site1
            ind2 <- pIDs != uID & lIDs == site2

            if (sum(ind1) != 0 & sum(ind2) != 0)
              if (dist.ave) {
                dist.o2 <- c(dist.o2, mean(dist.mat[ind1, ind2]))
              } else {
                dist.o2 <- c(dist.o2, dist.mat[ind1, ind2])
              }

          }

          if (!is.null(dist.o1) & !is.null(dist.o2)) {
            boxplot(list(ToSelf=dist.o1, ToOther=dist.o2), col='steelblue', ylab=paste(dist.name, 'distance'))
            stat.o <- mean(dist.o2) - mean(dist.o1)

            stat.p <- NULL
            for (i in 1:nPerm) {
              # Exchange the uterus between people :D
              pIDs.p <- pIDs
              pIDs.uterus <- pIDs[lIDs == site1]
              temp <- factor(pIDs.uterus)
              levels(temp) <- sample(levels(temp))
              pIDs.uterus.p <- as.character(temp)
              pIDs.p[lIDs == site1] <- pIDs.uterus.p

              # To increase the number of permutation
              pIDs.uterus <- pIDs[lIDs == site2]
              temp <- factor(pIDs.uterus)
              levels(temp) <- sample(levels(temp))
              pIDs.uterus.p <- as.character(temp)
              pIDs.p[lIDs == site2] <- pIDs.uterus.p

              # To self
              dist.o1 <- NULL
              for (uID in uIDs) {
                ind1 <- pIDs.p == uID & lIDs == site1
                ind2 <- pIDs.p == uID & lIDs == site2

                if (sum(ind1) != 0 & sum(ind2) != 0)
                  if (dist.ave) {
                    dist.o1 <- c(dist.o1, mean(dist.mat[ind1, ind2]))
                  } else {
                    dist.o1 <- c(dist.o1, dist.mat[ind1, ind2])
                  }
              }

              # To other people
              dist.o2 <- NULL
              for (uID in uIDs) {
                ind1 <- pIDs.p == uID & lIDs == site1
                ind2 <- pIDs.p != uID & lIDs == site2

                if (sum(ind1) != 0 & sum(ind2) != 0)
                  if (dist.ave) {
                    dist.o2 <- c(dist.o2, mean(dist.mat[ind1, ind2]))
                  } else {
                    dist.o2 <- c(dist.o2, dist.mat[ind1, ind2])
                  }
              }

              # boxplot(list(ToSelf=dist.o1, ToOther=dist.o2))
              stat.p <- c(stat.p, mean(dist.o2) - mean(dist.o1))
            }

            pv <- mean(c(stat.p >= stat.o, TRUE))

            hist(stat.p, xlab=paste(dist.name, 'distance'), col='steelblue', xlim=c(min(stat.p) * 2, max(stat.p) * 2),
                 main='Average distances observed under permutation', sub=paste('P < ', pv))

            abline(v=stat.o, col='red')
            pmat[dist.name, site1, site2] <- pmat[dist.name, site2, site1] <- pv
          }

        }
        dev.off()
      }

    }
  }
  cat('Finished!\n')

  pmat.csv <- NULL
  for (dist.name in dist.names) {
    temp <- as.matrix(pmat[dist.name, ,  ])
    mode(temp) <- 'character'
    pmat.csv <- rbind(pmat.csv, rbind(Distance=c(dist.name, rep('', ncol(temp) - 1)), temp, ''))
  }

  write.csv(pmat.csv, paste0('SiteCorrelationP_', ann, '.csv'))
  return(invisible(pmat))
}

# New: 2017_11_03
perform_site_taxa_correlation_test <- function (data.obj, site.name, sites, subject,
                                                taxa.levels = c('Genus'), taxa.names = NULL, dist.func = function (x) {dist(sqrt(x))},
                                                prev=0.1, minp=0.002, medianp=NULL, ann = '', ...) {
  pmat.list <- list()
  for (taxa.level in taxa.levels) {
    cat(taxa.level, '...\n')
    genus <- data.obj$abund.list[[taxa.level]]
    genus <- t(t(genus) / colSums(genus))
    if (is.null(taxa.names)) {
      ind <- rep(TRUE, nrow(genus))
      if (!is.null(prev)) ind <- ind & rowMeans(genus != 0) > prev
      if (!is.null(minp)) ind <- ind & rowMaxs(genus) > minp
      if (!is.null(medianp)) {
        nz.mean <- apply(genus, 1, function(x) median(x[x!=0]))
        ind <- ind & nz.mean > medianp
      }
      genus <- genus[ind, , drop = FALSE]
      gen.names <- rownames(genus)
    } else {
      if (sum(!(taxa.names %in% rownames(genus))) != 0) {
        stop('Could not find the taxa! Please check!\n')
      }
      genus <- genus[taxa.names, , drop = FALSE]
      gen.names <- taxa.names
    }

    # Convert into distance
    dist.obj <- list()
    for (gen.name in gen.names) {
      dist.obj[[gen.name]] <- as.matrix(dist.func(genus[gen.name, ]))
    }

    pmat.list[[taxa.level]] <- perform_site_correlation_test(data.obj, dist.obj, site.name, sites, subject,
                                                             dist.names = gen.names,  ann = paste(taxa.level, ann, sep = '_'), ...)

  }
  return(invisible(pmat.list))
}

# Rev: 2017_02_19 Add hei and wid
# Rev: 2018_01_15 Add qt.outlier=0.097
# Rev: 2018_03_06 position_jitter(h = 0), coord_cartesian
# Rev: 2018_07_16
generate_taxa_boxplot <- function (data.obj,  grp.name, strata=NULL, scale='P',  taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
                                   taxa.name='All', pseudo.ct=0.5, rm.outlier=TRUE, qt.outlier=0.97, prev=0.1, minp=0.002, ann='All',
                                   subject=NULL, l.size=0.5, p.size=2.5, hei=NULL, wid=NULL) {
  # To be completed
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  if (!is.null(subject)) {
    ID <- df[, subject]
  }

  if (is.null(hei)) {
    hei <- 5
  } else {
    hei <- hei
  }

  if (is.null(strata)) {
    if (is.null(wid)) {
      wid <- 5
    } else {
      wid <- wid
    }
  } else {
    if (is.null(wid)) {
      wid <- 5.5
    } else {
      wid <- wid
    }
  }

  if (scale == 'P') {
    ylab <- 'Proportion'
  }
  if (scale == 'logP') {
    ylab <- 'log10(Proportion)'
  }
  if (scale == 'sqrtP') {
    ylab <- 'sqrt(Proportion)'
  }
  if (scale == 'binary') {
    ylab <- 'Count'
  }

  for (LOI in taxa.levels) {
    if (LOI == 'All') {
      if (taxa.name == 'All') stop('Taxa names are not required to be all when taxa level is also all!\n')
      headnames <- NULL
      prop <- NULL
      for (LOI2 in names(data.obj$abund.list)) {
        prop2 <- data.obj$abund.list[[LOI2]]
        taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
        if (length(taxa.name2) != 0) {
          if (scale == 'logP') {
            prop2 <- prop2 + pseudo.ct
          }
          prop2 <- t(t(prop2) / colSums(prop2))
          headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
          names(headnames2) <- rownames(prop2)
          prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
          headnames <- c(headnames, headnames2)
        }
      }

    } else {

      cat(LOI, "\n")
      ct <- data.obj$abund.list[[LOI]]
      prop <- t(t(ct) / colSums(ct))
      ind <- rowMaxs(prop) > minp & rowSums(prop != 0) > prev * ncol(prop)

      if (scale == 'logP') {
        ct <- ct + pseudo.ct
        prop <- t(t(ct) / colSums(ct))
      }

      headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
      names(headnames) <- rownames(prop)

      if (taxa.name == 'All') {
        prop <- prop[ind, , drop=FALSE]
      } else {
        # Rev: 2018_07_16
        prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
      }


    }

    if (scale == 'logP') {
      prop <- log10(prop)
    }

    if (scale == 'sqrtP') {
      prop <- sqrt(prop)
    }

    if (scale == 'binary') {
      temp <- prop != 0
      prop[temp] <- 'Presence'
      prop[!temp] <- 'Absence'
    }


    pdf(paste("Taxa_Boxplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)

    if (is.null(strata)) {
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        if (scale != 'binary') {
          if (scale == 'P') {
            if (rm.outlier == T) {
              ylims <- c(min(taxon.abund),  min(quantile(taxon.abund, qt.outlier) * 1.25, 1))
            } else {
              ylims <- range(taxon.abund)
            }

          } else {
            ylims <- range(taxon.abund)
          }

          if (is.null(subject)) {
            df2 <- data.frame(Value=taxon.abund, Group=grp)
            dodge <- position_dodge(width=0.9)
            obj <- ggplot(df2, aes(x=Group, y=Value, col=Group)) +
              geom_boxplot(position=dodge, outlier.colour = NA, alpha=0.75, lwd=0.35) +
              geom_jitter(alpha=0.6, size=p.size,  position = position_jitter(w = 0.1, h = 0)) +
              labs(y=ylab, title=headnames[taxon])
            if (rm.outlier) {
              obj <- obj + coord_cartesian(ylim = ylims)
            }
            #		ylim(ylims[1], ylims[2]) +
            obj <- obj + theme(legend.position="none")
            print(obj)
          } else {
            df2 <- data.frame(Value=taxon.abund, Group=grp, subject=ID)
            dodge <- position_dodge(width=0.9)
            obj <- ggplot(df2, aes(x=Group, y=Value, shape=Group, group=subject)) +
              geom_point(size=p.size) +
              geom_line(size=l.size) +
              labs(y=ylab, title=headnames[taxon]) +
              #	ylim(ylims) +
              theme(legend.position="none")
            print(obj)

          }

        } else {
          df2 <- data.frame(Value=taxon.abund, Group=grp)

          obj <- ggplot(df2, aes(x=Group, fill=Value)) +
            geom_bar(width=.5) +
            labs(y=ylab, title=headnames[taxon]) +
            theme(legend.title=element_blank())

          print(obj)
        }

      }
    } else {
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        if (scale != 'binary') {
          if (scale == 'P') {
            if (rm.outlier == T) {
              ylims <- c(min(taxon.abund),  min(quantile(taxon.abund, qt.outlier) * 1.25, 1))
            } else {
              ylims <- range(taxon.abund)
            }

          } else {
            ylims <- range(taxon.abund)
          }
          grp2 <- df[, strata]
          df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
          dodge <- position_dodge(width=0.75)
          obj <- ggplot(df2, aes(x=Strata, y=Value, col=Group, fill=Group)) +
            geom_boxplot(fill='white',  position=dodge, outlier.colour = NA, alpha=0.75, width=0.65, lwd=0.35) +
            geom_point(position=position_jitterdodge(dodge.width=0.75), size=p.size, alpha=0.6) +
            labs(y=ylab, x=strata, title=headnames[taxon])
          if (rm.outlier) {
            obj <- obj + coord_cartesian(ylim = ylims)
          }
          #		ylim(ylims[1], ylims[2]) +
          obj <- obj + theme(legend.title=element_blank())
          print(obj)
        } else {
          grp2 <- df[, strata]
          df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
          obj <- ggplot(df2, aes(x=Group, fill=Value)) +
            geom_bar(width=.5) +
            labs(y=ylab, title=headnames[taxon]) +
            facet_wrap(~ Strata) +
            theme(legend.title=element_blank())

          print(obj)
        }
      }

    }
    dev.off()
  }

}

# Rev: 2017_02_19 Add Pseudo.ct, hei0, wid0; strata2; smooth.method
# Rev: 2018_01_15 Add qt.outlier
# Rev: 2018_07_16 Add an option to combine strata in the same plot, coord_cartesian
# Rev: 2019_01_03 Selecting taxa before adding pseudo-ct
# Note: if rm.outlier=TRUE, the smooth function only acts on the points in the plot region. May not be desirable
generate_taxa_scatterplot <- function (data.obj,  grp.name, subject=NULL,
                                       taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'), taxa.name='All', prev=0.1, minp=0.002,
                                       strata=NULL, strata2=NULL, combine.strata=FALSE, combine.strata2=FALSE, smooth.method='loess',
                                       pt.shape=16, pt.alpha=0.5, pt.size=2, subject.pt=FALSE,
                                       scale='P', pseudo.ct=0.5, is.pvalue=FALSE, rm.outlier=FALSE, qt.outlier=0.97,
                                       ann='All', hei=NULL, wid=NULL) {
  # To be completed
  df <- data.obj$meta.dat
  grp <- df[, grp.name]

  if(is.null(subject)) {
    ID <- NA
  } else {
    ID <- factor(df[, subject])
  }

  if (is.null(hei) | is.null(wid)) {
    hei <- 5
    wid <- 6
    if (scale != 'binary') {
      if (!is.null(strata) & combine.strata == FALSE) {
        wid <- 4 * nlevels(factor(df[, strata]))
      }
      if (!is.null(strata2) & combine.strata2 == FALSE) {
        hei <- 2.5 * nlevels(factor(df[, strata2]))
      }
    } else {
      if (!is.null(strata)) {
        wid <- 2.5 * nlevels(factor(df[, strata]))
      }
      if (!is.null(strata2)) {
        hei <- 2.5 * nlevels(factor(df[, strata2]))
      }
    }

  }

  if (scale == 'P') {
    ylab <- 'Proportion'
  }
  if (scale == 'logP') {
    ylab <- 'log10(Proportion)'
  }
  if (scale == 'sqrtP') {
    ylab <- 'sqrt(Proportion)'
  }
  if (scale == 'binary') {
    ylab <- 'Count'
  }

  for (LOI in taxa.levels) {
    if (LOI == 'All') {
      if (taxa.name == 'All') stop('Taxa names need to be specified when taxa level is  all!\n')
      headnames <- NULL
      prop <- NULL
      for (LOI2 in names(data.obj$abund.list)) {
        prop2 <- data.obj$abund.list[[LOI2]]
        taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
        if (length(taxa.name2) != 0) {
          if (scale == 'logP') {
            prop2 <- prop2 + pseudo.ct
          }
          prop2 <- t(t(prop2) / colSums(prop2))
          headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
          names(headnames2) <- rownames(prop2)
          prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
          headnames <- c(headnames, headnames2)
        }
      }

    } else {

      cat(LOI, "\n")
      ct <- data.obj$abund.list[[LOI]]
      prop <- t(t(ct) / colSums(ct))
      ind <- rowMaxs(prop) > minp & rowSums(prop != 0) > prev * ncol(prop)

      if (scale == 'logP') {
        ct <- ct + pseudo.ct
        prop <- t(t(ct) / colSums(ct))
      }

      headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
      names(headnames) <- rownames(prop)

      if (taxa.name == 'All') {
        prop <- prop[ind, , drop=FALSE]
      } else {
        prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
      }

    }

    if (scale == 'logP') {
      prop <- log10(prop)
    }

    if (scale == 'sqrtP') {
      prop <- sqrt(prop)
    }
    if (scale == 'asinsqrtP') {
      prop <- asin(sqrt(prop))
    }

    if (scale == 'binary') {
      temp <- prop != 0
      prop[temp] <- 'Presence'
      prop[!temp] <- 'Absence'
    }


    pdf(paste("Taxa_Scatterplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
    if (is.null(strata)) {
      if (is.pvalue == TRUE) {
        if (scale != 'binary') {
          pvalues <- apply(prop, 1, function (x) {
            formatC(cor.test(as.numeric(x), grp, method='spearman')$p.value, digits=3)
          })
        } else {
          pvalues <- apply(prop, 1, function (x) {
            if (nlevels(factor(x)) > 1) {
              return(formatC(wilcox.test(grp ~ x)$p.value, digits=3))
            } else {
              return('1.0')
            }

          })
        }

      }
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        if (scale != 'binary') {
          if (scale == 'P') {
            if (rm.outlier == T) {
              ylims <- c(min(taxon.abund),  min(quantile(taxon.abund, qt.outlier) * 1.25, 1))
            } else {
              ylims <- range(taxon.abund)
            }
          } else {
            ylims <- range(taxon.abund)
          }
          df2 <- data.frame(Value=taxon.abund, Group=grp, ID=ID)
          dodge <- position_dodge(width=0.9)
          obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
            #		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
            geom_smooth(method=smooth.method, aes(group=NA), size=0.75)
          if (is.pvalue == TRUE) {
            obj <- obj + labs(y=ylab, x=grp.name, title=paste0(headnames[taxon], '(P=', pvalues[taxon], ')'))
          } else {
            obj <- obj + labs(y=ylab, x=grp.name, title=paste0(headnames[taxon]))
          }
          if (rm.outlier) {
            obj <- obj + coord_cartesian(ylim = ylims)
          }
          obj <- obj + theme(legend.position="none")
          if (!is.null(subject)) {
            if (subject.pt == FALSE) {
              obj <- obj + geom_line(size=0.25, alpha=0.75)
            }  else {
              obj <- obj + geom_line(size=0.25, alpha=0.75) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha)
            }

          } else {
            obj <- obj + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha)
          }
          print(obj)
        } else {
          df2 <- data.frame(Value=taxon.abund, Group=grp, ID=ID)
          dodge <-
            obj <- ggplot(df2, aes(y=Group, x=Value)) +
            geom_boxplot(position=position_dodge(width=0.9), outlier.colour = NA, alpha=0.75) +
            geom_jitter(alpha=pt.alpha, size=pt.size,  position = position_jitter(w = 0.1, h = 0))

          if (is.pvalue == TRUE) {
            obj <- obj + labs(y=grp.name, x='Presence/Absence', title=paste0(headnames[taxon], '(P=', pvalues[taxon], ')'))
          } else {
            obj <- obj + labs(y=grp.name, x='Presence/Absence', title=paste0(headnames[taxon]))
          }
          obj <- obj + theme(legend.position="none")


          print(obj)
        }

      }
    } else {
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        if (scale != 'binary') {
          if (scale == 'P') {
            if (rm.outlier == T) {
              ylims <- c(min(taxon.abund), min(quantile(taxon.abund, qt.outlier) * 1.25, 1))
            } else {
              ylims <- range(taxon.abund)
            }

          } else {
            ylims <- range(taxon.abund)
          }
          grp2 <- factor(df[, strata])
          if (!is.null(strata2)) {
            grp3 <- factor(df[, strata2])
          } else {
            grp3 <- grp2
          }
          grp4 <- factor(paste(grp2, grp3))

          # VERY bad names with only difference in capital letters, variable names are also confusing
          # Rev: 2018_07_16

          df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2, Strata2=grp3, Strata3=grp4, ID=ID)


          if (is.null(strata2)) {
            if (combine.strata == TRUE) {
              obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata)) +
                #			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
                geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
                labs(y=ylab, x=grp.name, title=headnames[taxon])
            } else {
              obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
                #			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
                geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
                labs(y=ylab, x=grp.name, title=headnames[taxon]) +
                facet_grid(. ~ Strata)
            }

          } else {
            if (combine.strata == TRUE) {
              if (combine.strata2 == TRUE) {
                obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata3)) +
                  #			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
                  geom_smooth(method=smooth.method, aes(group=Strata3), size=0.75) +
                  labs(y=ylab, x=grp.name, title=headnames[taxon])
              } else {
                obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata)) +
                  #			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
                  geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
                  labs(y=ylab, x=grp.name, title=headnames[taxon]) +
                  facet_grid(Strata2 ~ .)
              }
            } else {
              if (combine.strata2 == TRUE) {
                warning('Currently, combine.strata2 will be forced to be FALSE if combine.strata=FALSE!\n')
              }
              obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
                #		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
                geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
                labs(y=ylab, x=grp.name, title=headnames[taxon]) +
                facet_grid(Strata2 ~ Strata)
            }

          }

          if (rm.outlier) {
            obj <- obj + coord_cartesian(ylim = ylims)
          }

          obj <- obj + theme(legend.title=element_blank())

          if (!is.null(subject)) {
            if (subject.pt == FALSE) {
              obj <- obj + geom_line(size=0.25, alpha=0.75)
            }  else {
              obj <- obj + geom_line(size=0.25, alpha=0.75) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha)
            }

          } else {
            obj <- obj + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha)
          }
          print(obj)

        } else {
          grp2 <- df[, strata]
          if (!is.null(strata2)) {
            grp3 <- df[, strata2]
          } else {
            grp3 <- grp2
          }
          df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2, Strata2=grp3)
          dodge <- position_dodge(width=0.9)

          if (is.null(strata2)) {
            obj <- ggplot(df2, aes(x=Strata, y=Group, col=Value, fill=Value)) +
              geom_boxplot(fill='white', position=dodge, outlier.colour = NA, alpha=0.75) +
              geom_point(position=position_jitterdodge(dodge.width=0.9), size=pt.size, alpha=pt.alpha) +
              labs(y=grp.name, x=strata, title=headnames[taxon])
          }  else {
            obj <- ggplot(df2, aes(x=Strata, y=Group, col=Value, fill=Value)) +
              geom_boxplot(fill='white', position=dodge, outlier.colour = NA, alpha=0.75) +
              geom_point(position=position_jitterdodge(dodge.width=0.9), size=pt.size, alpha=pt.alpha) +
              facet_grid(Strata2 ~ .) +
              labs(y=grp.name, x=strata, title=headnames[taxon])
          }

          print(obj)
        }
      }

    }
    dev.off()
  }
}

# New: 2019_01_03 - e.g., study the correlation between different sampling methods
# Restrict to two-level factor and pairs - labelling outlier has not been completed
generate_taxa_correlationplot <- function(data.obj,  grp.name, subject.name = NULL, strata = NULL,
                                          scale='sqrt', pseudo.ct = 0.5, error = 'se', n.bt = 50,  type = c('logP', 'Cor'), p.cutoff = 0.01,
                                          taxa.levels=c('Species'), taxa.name = 'All', prev = 0.1, minp = 0.002,
                                          k.outlier = 5, err.wid = 0, err.hei = 0, xlab = 'Grp1.mean', ylab = 'Grp2.mean',
                                          wid = 7.5, hei = 6, pt.size = 2,   ann = '') {

  type <- match.arg(type)
  # Strata = to color strata
  # We will only focus on pairs
  if (type == 'Cor' & is.null(subject.name)) stop('Correlation measure requires paired samples!\n')
  df <- data.obj$meta.dat
  grp <- factor(df[, grp.name])
  grp2 <- as.numeric(grp)
  ord <- 1:length(grp)

  if (!is.null(subject.name)) {
    subject <- factor(df[, subject.name])

    # Get one sample for each subject in each group
    ind <- as.vector(tapply(rownames(df), list(grp, subject), function (x) {
      if (length(x) < 1) {
        return(NULL)
      } else {
        sample(x, 1)
      }
    }))
    ind <- ind[!is.na(ind)]
    data.obj <- subset_data(data.obj, ind)

    df <- data.obj$meta.dat
    grp <- factor(df[, grp.name])
    subject <- factor(df[, subject.name])

    tab <- table(subject)
    ind <- subject %in% names(tab)[tab == 2]

    data.obj <- subset_data(data.obj, ind)
    df <- data.obj$meta.dat
    grp <- factor(df[, grp.name])
    subject <- factor(df[, subject.name])

    ord <- order(subject)
    df <- df[ord, ]
    grp <- grp[ord]
    grp2 <- as.numeric(grp)
    subject <- subject[ord]
  }

  plot.data.list <- list()
  for (i in 1:length(taxa.levels)) {
    LOI <- taxa.levels[i]
    cat(LOI, "\n")

    ct <- data.obj$abund.list[[LOI]][, ord]

    prop <- t(t(ct) / colSums(ct))
    ind <- rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop)

    if (scale == 'log') {
      ct <- ct + pseudo.ct
      prop <- t(t(ct / colSums(ct)))
    }

    if (taxa.name == 'All') {
      prop <- prop[ind, , drop=FALSE]
    } else {
      # Rev: 2018_07_16
      prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
    }

    if (is.null(strata)) {

      # Calculate spearman correlation
      df2 <- t(apply(prop, 1, function (x) {
        x1 <- x[grp2 == 1]
        x2 <- x[grp2 == 2]
        if (error == 'se') scale.factor <- 1
        if (error == 'ci') scale.factor <- 1.96
        se.x1 <- sd(sapply(1:n.bt, function (i) mean(sample(x1, repl=TRUE)))) * scale.factor
        se.x2 <- sd(sapply(1:n.bt, function (i) mean(sample(x2, repl=TRUE)))) * scale.factor
        m.x1 <- mean(x1)
        m.x2 <- mean(x2)
        if (!is.null(subject.name)) {
          cor.val <- cor(x1, x2, method = 'spearman')
          p.val <- wilcox.test(x1, x2, pair = TRUE)$p.value
        } else {
          cor.val <- NA
          p.val <- wilcox.test(x1, x2)$p.value
        }
        xmin <- m.x1 - se.x1
        xmax <- m.x1 + se.x1
        ymin <- m.x2 - se.x2
        ymax <- m.x2 + se.x2
        if (xmin <= 0) xmin <- 0.00001
        if (xmax >= 1) xmax <- 1 - 0.00001
        if (ymin <= 0) ymin <- 0.00001
        if (ymax >= 1) ymax <- 1 - 0.00001
        p.val.b <- as.numeric(p.val <= p.cutoff)

        c(m.x1, m.x2, xmin, xmax, ymin, ymax, cor.val, log10(p.val), p.val.b)
      }))

      df2 <- data.frame(rownames(prop), df2)
      colnames(df2) <- c('Taxa', 'x', 'y', 'xmin', 'xmax', 'ymin', 'ymax', 'Cor', 'Pval', 'Pval.b')
      df2[, 'Pval.b'] <- factor(ifelse(df2[, 'Pval.b'] == 1, paste0('p <=', p.cutoff), paste0('p >', p.cutoff)),
                                levels = c(paste0('p >', p.cutoff), paste0('p <=', p.cutoff)))
      #df2[, 'Pval.b'] <- factor(df2[, 'Pval.b'] )
      #df2 <- df2[!is.na(df2$Cor), ]

      obj1 <-  ggplot(data = df2, aes(x = x,y = y)) +
        geom_errorbar(aes(ymin = ymin, ymax = ymax), width = err.wid, color = 'gray') +
        geom_errorbarh(aes(xmin = xmin, xmax = xmax), height = err.hei, color = 'gray')
      if (type == "Cor" & !is.null(subject.name)) {
        obj1 <- obj1 + geom_point(aes(fill = Cor, shape = Pval.b), size = pt.size)
      }
      if (type == "logP") {
        obj1 <- obj1 + geom_point(aes(fill = Pval, shape = Pval.b), size = pt.size)
      }
      obj1 <- obj1 +
        geom_abline(intercept=0, slope=1, color = 'gray') +
        scale_shape_manual(values = c(21, 24)) +
        xlab(xlab) +
        ylab(ylab) +
        labs(fill= ifelse(type == 'logP', 'log10(P-val)', 'Correlation'), shape = NULL)

      if (scale == 'sqrt') {
        # obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
        obj1 <- obj1 + scale_y_sqrt(
          breaks = scales::trans_breaks("sqrt", function(x) x^2),
          labels = scales::trans_format("sqrt", scales::math_format(.x^2))) +
          scale_x_sqrt(
            breaks = scales::trans_breaks("sqrt", function(x) x^2),
            labels = scales::trans_format("sqrt", scales::math_format(.x^2)))
      }

      # To be revised
      if (scale == 'log') {
        obj1 <- obj1 + scale_y_log10(
          breaks = scales::trans_breaks("log10", function(x) 10^x),
          labels = scales::trans_format("log10", scales::math_format(10^.x))) +
          scale_x_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x)))
      }
      if (scale == 'boxcox') {
        obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3)) +
          scale_x_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
      }

    }

    plot.data.list[[LOI]] <- df2
    pdf(paste("Taxa_Correlation_", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
    print(obj1)
    dev.off()
  }
  return(invisible(plot.data.list))
}

# Rev: 2018_03_08 add bootstrap standard error (normal approximation)
# Rev: 2019_01_03 handling zeros for log scales
taxa_barplot_aggregate <- function (prop, df, grp.name, strata=NULL, scale='sqrt',  xsize=10, ylab='Proportion', error='ci') {

  if (scale == 'log') {
    # Half of the non-zero proportion
    prop <- t(apply(prop, 1, function (x) {
      temp <- min(x[x != 0])
      x[x == 0] <- temp / 2
      x
    }))
  }

  grp <- factor(df[, grp.name])

  if (is.null(strata)) {
    df2 <- data.frame(Group=grp, t(prop))
    df2 <- melt(df2)
    colnames(df2) <- c('Group', 'Taxa', 'Value')

    # Could be revised
    temp1 <- aggregate(Value ~ Group + Taxa, df2, mean)

    if (error == 'se') {
      temp2 <- aggregate(Value ~ Group + Taxa, df2, function(x) {
        #sd(x) / sqrt(length(x))
        sd(sapply(1:50, function(i) mean(sample(x, repl=TRUE))))
      })
    }
    if (error == 'ci') {
      temp2 <- aggregate(Value ~ Group + Taxa, df2, function(x) {
        #sd(x) / sqrt(length(x))
        2 * sd(sapply(1:50, function (i) mean(sample(x, repl=TRUE))))
      })
    }


    df2 <- cbind(temp1, temp2[, 3])
    colnames(df2) <- c('Group', 'Taxa', 'Mean', 'SE')


    limits <- aes(ymax = Mean + SE, ymin = ifelse(Mean - SE > 0, Mean - SE, 0))
    dodge <- position_dodge(width=0.90)

    obj1 <- ggplot(df2, aes(x=Taxa, y=Mean, fill=Group)) +
      geom_bar(position=dodge, stat="identity", alpha=0.75) +
      geom_bar(position=dodge, stat="identity", alpha=0.75, colour="black", show.legend=FALSE, size=0.25) +
      geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
      labs(y=paste(ylab), x='') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
      theme(legend.position="top", legend.title=element_blank())

    if (scale == 'sqrt') {
      # obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
      obj1 <- obj1 + scale_y_sqrt(
        breaks = scales::trans_breaks("sqrt", function(x) x^2),
        labels = scales::trans_format("sqrt", scales::math_format(.x^2)))
    }

    # To be revised
    if (scale == 'log') {
      obj1 <- obj1 + scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x)))
    }
    if (scale == 'boxcox') {
      obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
    }

  } else {
    grp2 <- factor(paste(strata, df[, strata]), levels=paste(strata, levels(df[, strata])))
    df2 <- data.frame(Group=grp, Strata=grp2, t(prop))
    df2 <- melt(df2)
    colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Value')

    # Could be revised
    temp1 <- aggregate(Value ~ Group + Strata + Taxa, df2, mean)

    if (error == 'se') {
      temp2 <- aggregate(Value ~ Group + Strata +Taxa, df2, function(x) {
        #sd(x) / sqrt(length(x))
        sd(sapply(1:50, function(i) mean(sample(x, repl=TRUE))))
      })
    }
    if (error == 'ci') {
      temp2 <- aggregate(Value ~ Group + Strata + Taxa, df2, function(x) {
        #sd(x) / sqrt(length(x))
        2 * sd(sapply(1:50, function (i) mean(sample(x, repl=TRUE))))
      })
    }


    df2 <- cbind(temp1, temp2[, 4])
    colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Mean', 'SE')

    limits <- aes(ymax = Mean + SE, ymin = ifelse(Mean - SE > 0, Mean - SE, 0))
    dodge <- position_dodge(width=0.90)

    obj1 <- ggplot(df2, aes(x=Taxa, y=Mean, fill=Group)) +
      geom_bar(position=dodge, stat="identity", alpha=0.75) +
      geom_bar(position=dodge, stat="identity", alpha=0.75, colour="black", show.legend=FALSE, size=0.25) +
      geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
      facet_wrap(~Strata, ncol=1) +
      labs(y=paste(ylab), x='') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
      theme(legend.position="top", legend.title=element_blank())
    if (scale == 'sqrt') {
      # obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
      obj1 <- obj1 + scale_y_sqrt(
        breaks = scales::trans_breaks("sqrt", function(x) x^2),
        labels = scales::trans_format("sqrt", scales::math_format(.x^2)))
    }
    if (scale == 'log') {
      obj1 <- obj1 + scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x)))
    }
    if (scale == 'boxcox') {
      obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
    }

  }
  return(obj1)

}

# Rev: 2018_03_08 Add geom_point and jsize
# Rev: 2019_01_03 handling zeros for log scales
taxa_boxplot_aggregate <- function (prop, df, grp.name, strata=NULL, scale='none', xsize=10, jsize=NULL, ylab='Proportion') {


  grp <- factor(df[, grp.name])

  if (scale == 'log') {
    # Half of the non-zero proportion
    prop <- t(apply(prop, 1, function (x) {
      temp <- min(x[x != 0])
      x[x == 0] <- temp / 2
      x
    }))
  }

  if (is.null(jsize)) {
    nT <- nrow(prop)
    jsize <- 2 / (nT %/% 20 + 1)
  }

  if (is.null(strata)) {
    df2 <- data.frame(Group=grp, t(prop))
    df2 <- melt(df2)
    colnames(df2) <- c('Group', 'Taxa', 'Value')

    dodge <- position_dodge(width=0.88)

    obj1 <- ggplot(df2, aes(x=Taxa, y=Value, fill=Group)) +
      geom_boxplot(position=dodge, alpha = ifelse(jsize == 0, 0.75, 0.75), outlier.alpha=ifelse(jsize == 0, 0.5, 0),  lwd=0.25, fatten=1)
    if (jsize != 0) {
      obj1 <- obj1 + geom_point(position=position_jitterdodge(dodge.width=0.88), size=jsize, alpha=0.3)
    }
    obj1 <- obj1 +
      labs(y=paste(ylab), x='') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
      theme(legend.position="top", legend.title=element_blank())
    if (scale == 'sqrt') {
      # obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
      obj1 <- obj1 + scale_y_sqrt(
        breaks = scales::trans_breaks("sqrt", function(x) x^2),
        labels = scales::trans_format("sqrt", scales::math_format(.x^2)))
    }
    if (scale == 'log') {
      obj1 <- obj1 + scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x)))
    }
    if (scale == 'boxcox') {
      obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
    }

  } else {
    grp2 <- factor(paste(strata, df[, strata]), levels=paste(strata, levels(df[, strata])))
    df2 <- data.frame(Group=grp, Strata=grp2, t(prop))
    df2 <- melt(df2)
    colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Value')

    dodge <- position_dodge(width=0.88)

    obj1 <- ggplot(df2, aes(x=Taxa, y=Value, fill=Group)) +
      geom_boxplot(position=dodge, alpha = ifelse(jsize == 0, 0.75, 0.75),  outlier.alpha=ifelse(jsize == 0, 0.5, 0),  lwd=0.25, fatten=1)
    if (jsize != 0) {
      obj1 <- obj1 + geom_point(position=position_jitterdodge(dodge.width=0.88), size=jsize, alpha=0.2)
    }
    obj1 <- obj1 +
      #	geom_point(position=position_jitterdodge(dodge.width=0.95), size=jsize, alpha=0.3) +
      facet_wrap(~Strata, ncol=1) +
      labs(y=paste(ylab), x='') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
      theme(legend.position="top", legend.title=element_blank())
    if (scale == 'sqrt') {
      # obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
      obj1 <- obj1 + scale_y_sqrt(
        breaks = scales::trans_breaks("sqrt", function(x) x^2),
        labels = scales::trans_format("sqrt", scales::math_format(.x^2)))
    }
    # To be revised
    if (scale == 'log') {
      obj1 <- obj1 + scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x)))
    }
    if (scale == 'boxcox') {
      obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
    }

  }
  return(obj1)

}

gmd_biplot <- function(data.obj, dist.obj, grp.name, components=4){
  results = list()
  library(GMDecomp)
  tab <- data.obj$abund.list[["Genus"]]
  #tab <- data.obj$otu.tab
  delta = dist.obj$WUniFrac ^2
  N = dim(tab)[2]
  P = dim(tab)[1]
  one_N = t(t(rep(1, N)))
  J = diag(1, N) - (1/N) * one_N %*% t(one_N)
  H = -(1/2) * J %*% delta %*% J

  test_eig <- eigen(H)
  pos_eig <- which(test_eig$values>=0)
  eig_vals <- diag(test_eig$values[pos_eig])
  eig_vec <- test_eig$vectors[,pos_eig]

  H_star = eig_vec %*% eig_vals %*% t(eig_vec)

  test.gmd <- GMD(X=as.matrix(clr(t(prop.table(tab+1,2)))),
                  H=H_star,
                  Q=diag(1, P),
                  K=components)
  gmd.order = order(rowSums(test.gmd$V[,1:2]^2), decreasing=T)
  plot.index = gmd.order[1:3]

  results$screeplot <- screeplot(test.gmd)

  results$biplot <- biplot(fit=test.gmd, index=plot.index, names=plot.names,
                           sample.col=data.obj$meta.dat[,grp.name],
                           arrow.col='grey50')
  return(results)

}


generate_taxa_biplot <- function (data.obj, taxa, trans='sqrt', grp.name, ann='', ...) {

  grp <- data.obj$meta.dat[, grp.name]

  prop <- NULL
  for (LOI2 in names(data.obj$abund.list)) {
    ct <- data.obj$abund.list[[LOI2]]
    ct <- ct + 1
    prop0 <- t(t(ct) / colSums(ct))
    prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa), , drop=FALSE])
  }
  colnames(prop) <- colnames(prop0)
  if (nrow(prop) != length(taxa)) {
    warnings('Some taxa not found in abundance lists! Please check the names!\n')
  }

  if (trans == 'normal') 	prop <- t(apply(prop, 1, function(x) qqnorm(x, plot=F)$x))
  if (trans == 'log') prop <- log(prop)
  if (trans == 'sqrt') prop <- sqrt(prop)
  if (trans == 'rank') prop <- t(apply(prop, 1, rank))

  wine.pca <- prcomp(t(prop), scale. = TRUE)
  g <- ggbiplot::ggbiplot(wine.pca, obs.scale = 1, var.scale = 1,
                groups = grp, ellipse = TRUE, circle = FALSE, ...)
  g <- g + scale_color_discrete(name = '') + theme_bw()
  g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
  pdf(paste0("Taxa_Biplot_", ann, ".pdf"), height=6, width=6)
  print(g)
  dev.off()
}

# Rev: 2017_05_19: new formula for width
#      Add presence and absence bar
# Rev: 2018_07_16
generate_taxa_barplot_aggregate <- function (data.obj,  grp.name, strata=NULL, scale='sqrt', pseudo.ct = 0.5,
                                             taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
                                             taxa.name='All',  prev=0.1, minp=0.002, ann='All', wids=NULL, heis=NULL, xsize=c(14, 12, 10, 9, 8)) {
  # To be completed
  df <- data.obj$meta.dat
  grp <- df[, grp.name]

  for (i in 1:length(taxa.levels)) {
    LOI <- taxa.levels[i]
    cat(LOI, "\n")

    ct <- data.obj$abund.list[[LOI]]
    prop <- t(t(ct / colSums(ct)))
    ind <- rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop)

    if (scale == 'log') {
      ct <- ct + pseudo.ct
      prop <- t(t(ct) / colSums(ct))
    }

    if (taxa.name == 'All') {
      prop <- prop[ind, , drop=FALSE]
    } else {
      # Rev: 2018_07_16
      prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
    }

    # decreasing
    prop <- prop[rev(order(rowMeans(prop))), ]

    if (is.null(wids) | is.null(heis)) {
      wid <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
      wid <- sqrt(nlevels(grp) / 2) * wid
      hei <- 7
    } else {
      wid <- wids[i]
      hei <- heis[i]
    }

    pdf(paste("Taxa_Barplot_Aggregate", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
    obj1 <- taxa_barplot_aggregate(prop, df, grp.name, strata, scale, xsize[i])
    print(obj1)
    dev.off()
  }
}

# Rev: 2017_05_19: new formula for width
# Rev: 2018_03_08  jsize
# Rev: 2018_07_16
generate_taxa_boxplot_aggregate <- function (data.obj,  grp.name, strata=NULL, scale='sqrt', pseudo.ct = 0.5,
                                             taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
                                             taxa.name='All',  prev=0.1, minp=0.002, ann='All', wids=NULL, heis=NULL, xsize=c(14, 12, 10, 9, 8, 6, 4), jsize = 0) {
  # To be completed
  df <- data.obj$meta.dat
  grp <- factor(df[, grp.name])

  for (i in 1:length(taxa.levels)) {
    LOI <- taxa.levels[i]
    cat(LOI, "\n")

    ct <- data.obj$abund.list[[LOI]]
    prop <- t(t(ct) / colSums(ct))
    ind <- rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop)

    if (scale == 'log') {
      ct <- ct + pseudo.ct
      prop <- t(t(ct / colSums(ct)))
    }

    if (taxa.name == 'All') {
      prop <- prop[ind, , drop=FALSE]
    } else {
      # Rev: 2018_07_16
      prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
    }

    # decreasing
    prop <- prop[rev(order(rowMeans(prop))), ]

    if (is.null(wids) | is.null(heis)) {
      wid <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
      wid <- sqrt(nlevels(grp) / 2) * wid

      hei <- 7
    } else {
      wid <- wids[i]
      hei <- heis[i]
    }

    pdf(paste("Taxa_Boxplot_Aggregate", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
    obj1 <- taxa_boxplot_aggregate (prop, df, grp.name, strata, scale, xsize[i], jsize=jsize)
    print(obj1)
    dev.off()
  }
}

# Add combined barplot with error bar and presence/absence bar (currently presence/absence bar is in generate_taxa_boxplot
# Rev: 2018_01_15 Add qt.outlier=0.97
# Rev: 2018_03_06 coord_cartesian
# Rev: 2018_05_31 Add x axis label control
# Rev: 2018_07_16
# Rev: 2019_01_03 add sqrtP and selecting taxa before adding pseudo-count
generate_taxa_barplot <- function (data.obj,  grp.name, strata=NULL, scale='P', pseudo.ct = 0.5,
                                   taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
                                   taxa.name='All', rm.outlier=T, qt.outlier=0.97, prev=0.1, minp=0.002, ann='All', x.axis.lab=FALSE, x.axis.size = 6) {
  # To be completed
  df <- data.obj$meta.dat
  grp <- factor(df[, grp.name])

  for (LOI in taxa.levels) {
    if (LOI == 'All') {
      if (taxa.name == 'All') stop('Taxa names are not allowed to be all when taxa level is also all!\n')
      headnames <- NULL
      prop <- NULL
      for (LOI2 in names(data.obj$abund.list)) {
        prop2 <- data.obj$abund.list[[LOI2]]
        taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
        if (length(taxa.name2) != 0) {
          if (scale == 'logP') {
            prop2 <- prop2 + pseudo.ct
          }
          prop2 <- t(t(prop2) / colSums(prop2))
          headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
          names(headnames2) <- rownames(prop2)
          prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
          headnames <- c(headnames, headnames2)
        }
      }

    } else {
      cat(LOI, "\n")
      ct <- data.obj$abund.list[[LOI]]
      prop <- t(t(ct) / colSums(ct))
      ind <- rowMaxs(prop) > minp & rowSums(prop != 0) > prev * ncol(prop)

      if (scale == 'logP') {
        ct <- ct + pseudo.ct
        prop <- t(t(ct) / colSums(ct))
      }

      headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
      names(headnames) <- rownames(prop)

      if (taxa.name == 'All') {
        prop <- prop[ind, , drop=FALSE]
      } else {
        # Rev: 2018_07_16
        prop <- prop[rownames(prop) %in% taxa.name, , drop=FALSE]
      }

    }

    if (scale == 'logP') {
      prop <- log10(prop)
    }

    if (scale == 'sqrtP') {
      prop <- sqrt(prop)
    }


    hei <- 5
    if (is.null(strata)) {
      wid <- 5
    } else {
      wid <- 4 * nlevels(df[, strata])
    }

    if (scale == 'P') {
      ylab <- 'Proportion'
    }

    if (scale == 'logP') {
      ylab <- 'log10(Proportion)'
    }

    if (scale == 'sqrtP') {
      ylab <- 'sqrt(Proportion)'
    }

    pdf(paste("Taxa_Barplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
    if (is.null(strata)) {
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        taxon.abund2 <- taxon.abund
        if (scale == 'P') {
          if (rm.outlier == T) {
            ylims <- c(0, quantile(taxon.abund, qt.outlier) * 1.25)
            taxon.abund[taxon.abund > quantile(taxon.abund, qt.outlier) * 1.25] <- quantile(taxon.abund, qt.outlier) * 1.25
          } else {
            ylims <- c(0, max(taxon.abund))
          }

        } else {
          ylims <- range(taxon.abund)
        }
        # df2 <- data.frame(x=factor(1:length(taxon.abund)), Value=taxon.abund, Group=grp)
        df2 <- data.frame(x=colnames(prop), Value=taxon.abund, Group=grp)
        rownames(df2) <- colnames(prop)

        dodge <- position_dodge(width=0.99)
        obj <- ggplot(df2, aes(x=x, y=Value, fill=Group)) +
          geom_bar(position=dodge, stat='identity', alpha=0.75) +
          facet_grid(. ~ Group, scales='free_x', space="free") +
          labs(y=ylab, title=headnames[taxon])

        if (rm.outlier) {
          obj <- obj + coord_cartesian(ylim = ylims)
        }

        if (!x.axis.lab) {
          obj <- obj + theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())
        } else {
          obj <- obj + theme(axis.text.x = element_text(size = x.axis.size, angle = 90, vjust = 0.5, hjust=1))
        }

        obj <- obj  +
          xlab('') +

          theme(legend.position="none")
        mean_df <- data.frame(Group=levels(grp), yint = tapply(taxon.abund2, grp, mean))
        obj <- obj + geom_hline(data = mean_df, aes(yintercept = yint))
        mean_df <- data.frame(Group=levels(grp), yint = tapply(taxon.abund2, grp, median))
        obj <- obj + geom_hline(data = mean_df, aes(yintercept = yint), linetype=2)

        print(obj)
      }
    } else {
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        taxon.abund2 <- taxon.abund
        if (scale == 'P') {
          if (rm.outlier == T) {
            ylims <- c(0, quantile(taxon.abund, qt.outlier) * 1.25)
            taxon.abund[taxon.abund > quantile(taxon.abund, qt.outlier) * 1.25] <- quantile(taxon.abund, qt.outlier) * 1.25
          } else {
            ylims <- c(0, max(taxon.abund))
          }

        } else {
          ylims <- range(taxon.abund)
        }

        grp2 <- df[, strata]
        obj.list <- list()
        for (level in levels(grp2)) {
          ind <- grp2 %in% level
          # df2 <- data.frame(x=factor(1:length(taxon.abund[ind])), Value=taxon.abund[ind], Group=grp[ind])
          df2 <- data.frame(x=colnames(prop)[ind], Value=taxon.abund[ind], Group=grp[ind])
          rownames(df2) <- colnames(prop)[ind]

          dodge <- position_dodge(width=0.99)
          obj0 <- ggplot(df2, aes(x=x, y=Value, fill=Group)) +
            geom_bar(position=dodge, stat='identity', alpha=0.75) +
            facet_grid(. ~ Group, scales='free_x', space="free") +
            labs(y=ylab, title=headnames[taxon])

          if (rm.outlier) {
            obj0 <- obj0 + coord_cartesian(ylim = ylims)
          }

          if (!x.axis.lab) {
            obj0 <- obj0 + theme(axis.ticks.x=element_blank(), axis.text.x = element_blank())
          } else {
            obj0 <- obj0 + theme(axis.text.x = element_text(size = x.axis.size, angle = 90, vjust = 0.5, hjust=1))
          }


          obj0 <- obj0  +
            xlab(paste(strata, level)) +
            theme(legend.position="none")
          mean_df <- data.frame(Group=levels(grp[ind]), yint = tapply(taxon.abund2[ind], grp[ind], mean))
          obj0 <- obj0 + geom_hline(data = mean_df, aes(yintercept = yint))
          mean_df <- data.frame(Group=levels(grp[ind]), yint = tapply(taxon.abund2[ind], grp[ind], median))
          obj0 <- obj0 + geom_hline(data = mean_df, aes(yintercept = yint), linetype=2)
          obj.list[[level]] <- obj0
        }
        multiplot(plotlist=obj.list, cols=nlevels(grp2))
      }

    }
    dev.off()
  }
}


plot_effect_size <- function (month,  value, pos.lab, neg.lab, ylab, hjust1=1.3, hjust2=-0.3, lab.size=3, xsize=10) {
  month <- factor(month, levels=month)
  dtm <- data.frame(month=month, value=value)
  dtm$colour <- factor(ifelse(dtm$value < 0, neg.lab, pos.lab), levels=c(pos.lab, neg.lab))
  dtm$hjust <- ifelse(dtm$value > 0, hjust1, hjust2)
  obj <- ggplot(dtm, aes(month, value, label = month, hjust = hjust)) +
    geom_text(aes(y = 0, colour = colour), size=lab.size) +
    geom_bar(stat = "identity", aes(fill = colour)) +
    theme(axis.text.x = element_text(size=xsize)) +
    ylim(c(-max(abs(value))*1.1, max(abs(value))*1.1)) +
    coord_flip() +
    scale_x_discrete(breaks = NULL) +
    labs(x = "", y = ylab) +
    theme(legend.position="top", legend.title=element_blank())
}


plot_effect_size2 <- function (fold.dat.plot1, ylabel='log(Fold change)', is.ln=TRUE, ord=TRUE) {

  if (is.ln) {
    fold.dat.plot1[, c('Estimate', 'LCI', 'UCI')] <- fold.dat.plot1[, c('Estimate', 'LCI', 'UCI')] * log2(exp(1))
  }
  Alphas <- seq(1, 99, 2) / 100
  #Multiplier <- qnorm(1 - Alphas / 2)
  Multiplier <- seq(1, 0.01, len=50)
  zzTransparency <<- 1/(length(Multiplier)/4)

  #	fold.dat.plot1 <- data.frame(IV=rownames(fold.dat.plot), Estimate=fold.dat.plot[, 1], LCI=fold.dat.plot[, 2], UCI=fold.dat.plot[, 3])
  fold.dat.plot1 <- data.frame(cbind(fold.dat.plot1, Scalar=rep(Multiplier, each = nrow(fold.dat.plot1))))
  #	fold.dat.plot1[,  -1] <- apply(fold.dat.plot1[, -1], 2, function(x){as.numeric(as.character(x))})
  fold.dat.plot1$Emphasis <- by(1 - seq(0, 1, length = length(Multiplier) + 1)[-(length(Multiplier) + 1)],
                                as.character(round(Multiplier, 5)), mean)[as.character(round(fold.dat.plot1$Scalar, 5))]

  fold.dat.plot1$IV <- factor(fold.dat.plot1$IV, unique(fold.dat.plot1$IV))

  if (ord) {
    fold.dat.plot1 <- fold.dat.plot1[order(as.character(fold.dat.plot1$IV)), ]
  }

  OutputPlot <- ggplot2::qplot(data = fold.dat.plot1, x = IV, y = Estimate,
                               ymin = Estimate - (Estimate -LCI)*Scalar, ymax = Estimate + (UCI - Estimate)*Scalar,
                               ylab = NULL, xlab = NULL, alpha = I(zzTransparency), colour = I(gray(0)), geom = "blank")

  OutputPlot <- OutputPlot + geom_hline(yintercept = 0, lwd = I(7/12), colour = I(hsv(0/12, 7/12, 7/12)), alpha = I(5/12))
  OutputPlot <- OutputPlot + geom_linerange(data = fold.dat.plot1, aes(size = 1/Emphasis), alpha = I(zzTransparency), colour = I(gray(0)))
  OutputPlot <- OutputPlot + scale_size_continuous() + guides(size=FALSE)
  #OutputPlot <- OutputPlot + facet_grid(~ ModelName)
  OutputPlot <- OutputPlot + coord_flip() + geom_point(aes(x = IV, y = Estimate), colour = I(gray(0))) + theme_bw() + ylab(ylabel)
}


# Rev: 2016_09_13 add glmmPQL-based overdipersed Poisson and binomial model
# glmmPQL default has overdispersion parameter. Regular Poisson and Binomial does not work.
perform_differential_analysis_para_single_RE <- function (taxon, ldep, grp.name, adj.name=NULL, subject, df, method='NB', LRT=FALSE) {
  # ldep: log depth (size factor)
  if (!is.null(adj.name) ) {
    if (sum(grepl(grp.name, c(adj.name)))) {
      stop('grp.name could not be part of adj.name, or there will be problem!\n')
    }
  }

  if (is.null(subject)) {
    stop('Random effects model require the specification of the subject parameter!\n')
  }

  df$ldep <- ldep
  df$taxon <- taxon
  if (LRT) warning('Currently, only Wald test is implemented!\n')


  if (method %in% c('NB', 'ZINB1', 'ZINB0', 'B0')) {
    if (is.null(adj.name)) {
      grp.name.adj.name.subject <- grp.name.adj.name.subject <- paste(grp.name,  '+ (1|', subject, ')')
    } else {

      grp.name.adj.name.subject <- paste(grp.name, '+', adj.name, '+ (1|', subject, ')')
    }

    if (method == 'NB') {
      m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df,
                        zeroInflation=FALSE, family='nbinom')
    }

    # 'B0' uses the glmmadmb, which will be compared to 'B' method of glmmPQL
    if (method == 'B0') {
      taxon2 <- as.numeric(taxon != 0)
      df$taxon2 <- taxon2
      m1.nb <- glmmadmb(as.formula(paste('taxon2 ~', grp.name.adj.name.subject, '+ ldep')), data = df,
                        zeroInflation=FALSE, family='binomial')
    }

    if (method == 'ZINB1') {
      m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df,
                        zeroInflation=TRUE, family='nbinom')
    }

    if (method == 'ZINB0') {
      m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df,
                        zeroInflation=TRUE, family='truncnbinom')
    }

    code <- list(m1.conv=m1.nb$convmsg)
    pv.nb <- wald.test(b = coef(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
    method <- paste(method, 'Wald')


    aic.nb <- -2 * m1.nb$loglik + 2 * m1.nb$npar
    coef.nb <- coef(m1.nb)
    ci.nb <- confint.default(m1.nb)
  }

  if (method %in% c('OP', 'B', 'QB')) {
    if (is.null(adj.name)) {
      grp.name.adj.name <- grp.name
    } else {
      grp.name.adj.name <- paste(grp.name, '+', adj.name)
    }

    subject <- paste('~ 1 |', subject)
    if (method == 'OP') {
      m1.nb <- glmmPQL(as.formula(paste('taxon ~', grp.name.adj.name, '+ offset(ldep)')), random = as.formula(subject), data = df,
                       verbose=F, family=quasipoisson())
    }

    # 'B0' uses the glmmadmb, which will be compared to 'B' method of glmmPQL
    if (method == 'B') {
      taxon2 <- as.numeric(taxon != 0)
      df$taxon2 <- taxon2
      m1.nb <- glmmPQL(as.formula(paste('taxon2 ~', grp.name.adj.name, '+ ldep')),  random = as.formula(subject), data = df,
                       verbose=F, family=binomial)
    }

    if (method == 'QB') {
      taxon2 <- as.numeric(taxon != 0)
      df$taxon2 <- taxon2
      m1.nb <- glmmPQL(as.formula(paste('taxon2 ~', grp.name.adj.name, '+ ldep')),  random = as.formula(subject), data = df,
                       verbose=F, family=quasibinomial())
    }

    pv.nb <- wald.test(b = fixed.effects(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
    method <- paste(method, 'Wald')
    code <- NA   # placeholder
    aic.nb <- NA # NA placeholder
    coef.nb <- fixed.effects(m1.nb)
    ci.nb <- confint.lme(m1.nb)

  }

  fc.nb <- coef.nb[grep(grp.name, names(coef.nb))]
  obj <- ci.nb[grep(grp.name, rownames(ci.nb)), ]

  if (is.vector(obj)) {
    fc.lc.nb <- obj[1]
    fc.uc.nb <- obj[2]
  } else {
    fc.lc.nb <- obj[, 1]
    fc.uc.nb <- obj[, 2]
    names(fc.lc.nb) <- paste(names(fc.lc.nb), '2.5%')
    names(fc.uc.nb) <- paste(names(fc.uc.nb), '97.5%')
  }

  return(list(method=method, pv=pv.nb, lfc=fc.nb, lfc.lci=fc.lc.nb, lfc.uci=fc.uc.nb, aic=aic.nb, code=code))
}



perform_differential_analysis_para_single_FE <- function (taxon.abund, ldep, grp.name, adj.name=NULL, subject=NULL, df, method='NB', LRT=FALSE) {
  # ldep: log depth (size factor)
  if (!is.null(adj.name)) {
    if (sum(grepl(grp.name, c(adj.name)))) {
      stop('grp.name could not be part of adj.name or subject, or there will be problem!\n')
    }
  }

  if (!is.null(subject)) {
    warnings('Fixed effects model will ignore the subject variable! Please use randome effects model!\n')
  }
  if (LRT & method == 'OP') warning('Overdispersed Poisson does not support LRT! Wald test used!\n')

  if (is.null(adj.name)) {
    grp.name.adj.name <- grp.name
  } else {
    grp.name.adj.name <- paste(grp.name, '+', adj.name)
  }
  if (method == 'NB') {
    m1.nb <- glm.nb(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep)')), data = df)
    if (LRT) {
      m0.nb <- update(m1.nb, as.formula(paste('. ~ . -',  grp.name)))
      code <- list(m1.conv=m1.nb$converged, m1.bound=m1.nb$boundary, m0.conv=m0.nb$converged, m0.bound=m0.nb$boundary)
      pv.nb <- anova(m1.nb, m0.nb)['Pr(Chi)'][2, ]
      method <- paste(method, 'LRT')
    } else {
      code <- list(m1.conv=m1.nb$converged, m1.bound=m1.nb$boundary)
      pv.nb <- wald.test(b = coef(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
      method <- paste(method, 'Wald')
    }

    aic.nb <- summary(m1.nb)$aic

    coef.nb <- coef(m1.nb)
    fc.nb <- coef.nb[grep(grp.name, names(coef.nb))]

    ci.nb <- confint.default(m1.nb)
    obj <- ci.nb[grep(grp.name, rownames(ci.nb)), ]

    if (is.vector(obj)) {
      fc.lc.nb <- obj[1]
      fc.uc.nb <- obj[2]
    } else {
      fc.lc.nb <- obj[, 1]
      fc.uc.nb <- obj[, 2]
      names(fc.lc.nb) <- paste(names(fc.lc.nb), '2.5%')
      names(fc.uc.nb) <- paste(names(fc.uc.nb), '97.5%')
    }
    return(list(method=method, pv=pv.nb, lfc=fc.nb, lfc.lci=fc.lc.nb, lfc.uci=fc.uc.nb, aic=aic.nb, code=code))
  }
  if (method == 'B') {
    taxon.abund2 <- as.numeric(taxon.abund != 0)
    m1.b <- glm(as.formula(paste('taxon.abund2 ~', grp.name.adj.name, '+ ldep')), data = df, family=binomial)
    if (LRT) {
      m0.b <- update(m1.b, as.formula(paste('. ~ . -',  grp.name)))
      code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary, m0.conv=m0.b$converged, m0.bound=m0.b$boundary)
      pv.b <- pchisq(2 * (logLik(m1.b) - logLik(m0.b)), df = df.residual(m0.b) - df.residual(m1.b), lower.tail=FALSE)
      method <- paste(method, 'LRT')
    } else {
      code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary)
      pv.b <- wald.test(b = coef(m1.b), Sigma = vcov(m1.b), Terms = grep(grp.name, names(coef(m1.b))))$result$chi2['P']
      method <- paste(method, 'Wald')
    }

    aic.b <- summary(m1.b)$aic
    coef.b <- coef(m1.b)
    fc.b <- coef.b[grep(grp.name, names(coef.b))]

    ci.b <- confint.default(m1.b)
    obj <- ci.b[grep(grp.name, rownames(ci.b)), ]

    if (is.vector(obj)) {
      fc.lc.b <- obj[1]
      fc.uc.b <- obj[2]
    } else {
      fc.lc.b <- obj[, 1]
      fc.uc.b <- obj[, 2]
      names(fc.lc.b) <- paste(names(fc.lc.b), '2.5%')
      names(fc.uc.b) <- paste(names(fc.uc.b), '97.5%')
    }
    return(list(method=method, pv=pv.b, lfc=fc.b, lfc.lci=fc.lc.b, lfc.uci=fc.uc.b, aic=aic.b, code=code))
  }

  # Rev: 2016_09_13 add 'QB', No likelihood ratio test
  if (method == 'QB') {
    taxon.abund2 <- as.numeric(taxon.abund != 0)
    m1.b <- glm(as.formula(paste('taxon.abund2 ~', grp.name.adj.name, '+ ldep')), data = df, family=quasibinomial)

    code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary)
    pv.b <- wald.test(b = coef(m1.b), Sigma = vcov(m1.b), Terms = grep(grp.name, names(coef(m1.b))))$result$chi2['P']
    method <- paste(method, 'Wald')


    aic.b <- summary(m1.b)$aic
    coef.b <- coef(m1.b)
    fc.b <- coef.b[grep(grp.name, names(coef.b))]

    ci.b <- confint.default(m1.b)
    obj <- ci.b[grep(grp.name, rownames(ci.b)), ]

    if (is.vector(obj)) {
      fc.lc.b <- obj[1]
      fc.uc.b <- obj[2]
    } else {
      fc.lc.b <- obj[, 1]
      fc.uc.b <- obj[, 2]
      names(fc.lc.b) <- paste(names(fc.lc.b), '2.5%')
      names(fc.uc.b) <- paste(names(fc.uc.b), '97.5%')
    }
    return(list(method=method, pv=pv.b, lfc=fc.b, lfc.lci=fc.lc.b, lfc.uci=fc.uc.b, aic=aic.b, code=code))
  }

  if (method == 'OP') {
    # No LRT
    m1.op <- glm(as.formula(paste('taxon.abund ~', grp.name.adj.name)), offset=ldep, data = df, family=quasipoisson)
    code <- list(m1.conv=m1.op$converged, m1.bound=m1.op$boundary)

    # pv.op <- pchisq(2 * (logLik(m1.op) - logLik(m0.op)), df = df.residual(m0.op) - df.residual(m1.op), lower.tail=FALSE) # LRT not applicable
    coef.op <- coef(m1.op)
    pv.op <- wald.test(b = coef.op, Sigma = vcov(m1.op), Terms = grep(grp.name, names(coef.op)))$result$chi2['P']
    method <- paste(method, 'Wald')
    fc.op <- coef.op[grep(grp.name, names(coef.op))]

    ci.op <- confint.default(m1.op)
    obj <- ci.op[grep(grp.name, rownames(ci.op)), ]

    if (is.vector(obj)) {
      fc.lc.op <- obj[1]
      fc.uc.op <- obj[2]
    } else {
      fc.lc.op <- obj[, 1]
      fc.uc.op <- obj[, 2]
      names(fc.lc.op) <- paste(names(fc.lc.op), '2.5%')
      names(fc.uc.op) <- paste(names(fc.uc.op), '97.5%')
    }
    return(list(method=method, pv=pv.op, lfc=fc.op, lfc.lci=fc.lc.op, lfc.uci=fc.uc.op, aic=NULL, code=code))
  }

  if (method == 'ZINB0') {
    m1.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep)')),
                        data = df, dist = "negbin", EM = TRUE)
    if (LRT) {
      if (is.null(adj.name)) {
        m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep)')),
                            data = df, dist = "negbin", EM = TRUE)
      } else {
        m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep)')),
                            data = df, dist = "negbin", EM = TRUE)
      }
      code <- list(m1.conv=m1.zinb$converged, m0.conv=m0.zinb$converged)
      # LRT
      pv.zinb <-  pchisq(2 * (logLik(m1.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m1.zinb), lower.tail=FALSE)
      method <- paste(method, 'LRT')
    } else {
      code <- list(m1.conv=m1.zinb$converged)
      pv.zinb <- wald.test(b = coef(m1.zinb), Sigma = vcov(m1.zinb), Terms = grep(grp.name, names(coef(m1.zinb))))$result$chi2['P']
      method <- paste(method, 'Wald')
    }

    aic.zinb <- -2 * logLik(m1.zinb) + 2 * (m1.zinb$n - m1.zinb$df.residual)

    coef.zinb <- coef(m1.zinb)
    fc.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]

    ci.zinb <- confint.default(m1.zinb)
    obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]

    if (is.vector(obj)) {
      fc.lc.zinb <- obj[1]
      fc.uc.zinb <- obj[2]
    } else {
      fc.lc.zinb <- obj[, 1]
      fc.uc.zinb <- obj[, 2]
      names(fc.lc.zinb) <- paste(names(fc.lc.zinb), '2.5%')
      names(fc.uc.zinb) <- paste(names(fc.uc.zinb), '97.5%')
    }
    return(list(method=method, pv=pv.zinb, lfc=fc.zinb, lfc.lci=fc.lc.zinb, lfc.uci=fc.uc.zinb, aic=aic.zinb, code=code))
  }


  if (method == 'ZINB1') {
    m1.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep) | ldep')),
                        data = df, dist = "negbin", EM = TRUE)
    if (LRT) {
      if (is.null(adj.name)) {
        m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep) | ldep')),
                            data = df, dist = "negbin", EM = TRUE)
      } else {
        m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep) | ldep')),
                            data = df, dist = "negbin", EM = TRUE)
      }
      code <- list(m1.conv=m1.zinb$converged, m0.conv=m0.zinb$converged)
      # LRT
      pv.zinb <-  pchisq(2 * (logLik(m1.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m1.zinb), lower.tail=FALSE)
      method <- paste(method, 'LRT')
    } else {
      code <- list(m1.conv=m1.zinb$converged)
      pv.zinb <- wald.test(b = coef(m1.zinb), Sigma = vcov(m1.zinb), Terms = grep(grp.name, names(coef(m1.zinb))))$result$chi2['P']
      method <- paste(method, 'Wald')
    }

    aic.zinb <- -2 * logLik(m1.zinb) + 2 * (m1.zinb$n - m1.zinb$df.residual)

    coef.zinb <- coef(m1.zinb)
    fc.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]

    ci.zinb <- confint.default(m1.zinb)
    obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]

    if (is.vector(obj)) {
      fc.lc.zinb <- obj[1]
      fc.uc.zinb <- obj[2]
    } else {
      fc.lc.zinb <- obj[, 1]
      fc.uc.zinb <- obj[, 2]
      names(fc.lc.zinb) <- paste(names(fc.lc.zinb), '2.5%')
      names(fc.uc.zinb) <- paste(names(fc.uc.zinb), '97.5%')
    }
    return(list(method=method, pv=pv.zinb, lfc=fc.zinb, lfc.lci=fc.lc.zinb, lfc.uci=fc.uc.zinb, aic=aic.zinb, code=code))
  }


  if (method == 'ZINB2') {
    m2.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep) |', grp.name.adj.name, '+ ldep')),
                        data = df, dist = "negbin", EM = TRUE)
    if (LRT) {
      if (is.null(adj.name)) {
        m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep) | ldep')),
                            data = df, dist = "negbin", EM = TRUE)
      } else {
        m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep) |', adj.name, ' + ldep')),
                            data = df, dist = "negbin", EM = TRUE)
      }

      code <- list(m1.conv=m2.zinb$converged, m0.conv=m0.zinb$converged)
      # LRT
      pv2.zinb <-  pchisq(2 * (logLik(m2.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m2.zinb), lower.tail=FALSE)
      method <- paste(method, 'LRT')
    } else {
      code <- list(m2.conv=m2.zinb$converged)
      pv2.zinb <- wald.test(b = coef(m2.zinb), Sigma = vcov(m2.zinb), Terms = grep(grp.name, names(coef(m2.zinb))))$result$chi2['P']
      method <- paste(method, 'Wald')
    }

    aic2.zinb <- -2 * logLik(m2.zinb) + 2 * (m2.zinb$n - m2.zinb$df.residual)

    coef.zinb <- coef(m2.zinb)
    fc2.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]

    ci.zinb <- confint.default(m2.zinb)
    obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]

    if (is.vector(obj)) {
      fc2.lc.zinb <- obj[1]
      fc2.uc.zinb <- obj[2]
    } else {
      fc2.lc.zinb <- obj[, 1]
      fc2.uc.zinb <- obj[, 2]
      names(fc2.lc.zinb) <- paste(names(fc2.lc.zinb), '2.5%')
      names(fc2.uc.zinb) <- paste(names(fc2.uc.zinb), '97.5%')
    }
    return(list(method=method, pv=pv2.zinb, lfc=fc2.zinb, lfc.lci=fc2.lc.zinb, lfc.uci=fc2.uc.zinb, aic=aic2.zinb, code=code))
  }

}

# Rev: 2016_09_14 Add Adaptive 0: switch OP and QB depending on the number of zeros.
# Rev: 2017_10_30  Support filtering based on coefficient of variation
# Rev: 2018_01_30 Add ct.min and handle NA size factor
perform_differential_analysis_para <- function (data.obj,  grp.name, adj.name=NULL, subject=NULL, RE=FALSE, method='Adaptive0', zerop.cutoff=0.25, ZINB='ZINB1', LRT=FALSE,
                                                taxa.levels=c('Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'), winsor=TRUE, winsor.qt=0.97, norm='GMPR', norm.level='Species', intersect.no=4, ct.min = 2,
                                                prev=0.1, minp=0.002, medianp=NULL, cv=NULL, mt.method='fdr', cutoff=0.15, ann='', ...) {
  # To be completed
  # subject holds the random effects formula
  if (!RE) {
    if (!(method %in% c('ZINB', 'B', 'QB', 'NB', 'OP', 'Adaptive0', 'Adaptive1', 'Adaptive2'))) stop('The speficied model is not supported!\n')
    perform_differential_analysis_para_single <- perform_differential_analysis_para_single_FE
    if (!is.null(subject)) warning('subject will not be used. Are you sure you want to run fixed effects model? ')
  } else {
    if (!(method %in% c('ZINB', 'B', 'B0', 'QB', 'NB', 'OP', 'Adaptive0', 'Adaptive1', 'Adaptive2'))) stop('The speficied model does not have random effects implementation!\n')
    if (ZINB != 'ZINB1') stop('Currently only ZINB1 is supported!\n')
    if (is.null(subject)) warning('subject is not supplied. Fixed effects model will be used instead!\n')
    perform_differential_analysis_para_single <- perform_differential_analysis_para_single_RE
  }

  df <- data.obj$meta.dat
  grp <- df[, grp.name]

  ind <- !is.na(grp)
  data.obj <- subset_data(data.obj, ind)
  grp <- grp[ind]
  df <- df[ind, ]

  if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
    data.obj$abund.list[['Species']] <- data.obj$otu.tab
    rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":",
                                                         data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
  }

  dep <- colSums(data.obj$otu.tab)
  diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']

  if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
    cat("Signficant sequencing depth confounding!\n")
    cat("For parametric test with sequence depth adjustment, please be cautious about the results!\n")
    cat("There may be potential residual sequence depth confounding!\n")
  }

  pv.list <- qv.list <-  fc.list <- fc.lc.list <- fc.uc.list <- met.list <- list()
  res.final <- NULL

  if (norm == 'Precalculated') {
    dep <- data.obj$size.factor
  }

  if (norm == 'GMPR') {
    if (norm.level %in% c('OTU', 'Species')) {
      tab <- data.obj$otu.tab
    } else {
      tab <- data.obj$abund.list[[norm.level]]
    }
    dep <- GMPR(tab, intersect.no, ct.min)

    # Rev: 2018_01_30
    ind <- !is.na(dep)
    data.obj <- subset_data(data.obj, ind)
    grp <- grp[ind]
    df <- df[ind, ]
    dep <- dep[ind]

  }

  if (norm == 'TSS') {
    dep <- colSums(data.obj$otu.tab)
  }

  ldep <- log(dep)

  for (LOI in taxa.levels) {
    cat(LOI, "\n")

    taxon.ct <- data.obj$abund.list[[LOI]]


    if (winsor == TRUE) {
      # Addressing the outlier (97% percent) or at least one outlier

      taxon.ct.p <- t(t(taxon.ct) / dep)
      taxon.ct.p <- apply(taxon.ct.p, 1, function(x) {
        cutoff <- quantile(x, winsor.qt)
        x[x >= cutoff] <- cutoff
        x
      }
      )
      # column/row switch
      taxon.ct <- t(round(taxon.ct.p * dep))

    }

    prop <- t(t(taxon.ct) / colSums(taxon.ct))
    #		if (!is.null(minp)) {
    #			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
    #			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
    #		}
    #
    #		if (!is.null(medianp)) {
    #			nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
    #			prop <- prop[nz.mean > medianp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
    #			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
    #		}

    if (!is.null(prev)) {
      prop <- prop[rowSums(prop!=0) > prev * ncol(prop), , drop=FALSE]
      taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
    }
    if (!is.null(minp)) {
      prop <- prop[rowMaxs(prop) > minp, , drop=FALSE]
      taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
    }

    if (!is.null(medianp)) {
      nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
      prop <- prop[nz.mean > medianp, , drop=FALSE]
      taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
    }

    if (!is.null(cv)) {
      prop <- prop[rowSds(prop) / rowMeans(prop) > cv, , drop=FALSE]
      taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
    }

    pv.vec <-  fc.vec <- fc.lc.vec <- fc.uc.vec <- met.vec <- conv.vec <- NULL
    obj <- NULL
    for (taxon in rownames(taxon.ct)) {
      cat('.')
      taxon.abund <- taxon.ct[taxon, ]

      ######## Logistic regression ###############
      if (method == 'B0') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='B0', LRT, ...))
      if (method == 'B') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='B', LRT, ...))
      if (method == 'QB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='QB', LRT, ...))
      ######## Overdispersed Poisson regression #########
      if (method == 'OP') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='OP', LRT, ...))
      ######## Negative binomial regression #########
      if (method == 'NB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...))
      ######## Zeroinflated negbinomial regression 1 ########
      if (method == 'ZINB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...))

      # Adpative 0 selects OP and QB based on the zero proportion (Not optimal)
      if (method == 'Adaptive0') {
        temp <- mean(as.numeric(taxon.abund != 0))

        if (temp > zerop.cutoff) {
          error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='QB', LRT, ...))
        } else {
          error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='OP', LRT, ...))
        }
      }

      # Adpative 1 selects NB and ZIB based on AIC
      if (method == 'Adaptive1') {
        error1 <- try(obj1 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...))
        error2 <- try(obj2 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...))
        if (class(error1) != 'try-error' & class(error2) != 'try-error') {
          if (obj1$aic < obj2$aic) {
            obj <- obj1
          } else {
            obj <- obj2
          }
          error <- error1
        } else {
          # pv == 0 indicates some problems in fitting
          if (class(error1) != 'try-error' & obj1$pv != 0) {
            obj <- obj1
            error <- error1
          } else {
            if (class(error2) != 'try-error' & obj1$pv != 0) {
              obj <- obj2
              error <- error2
            } else {
              error <- error2
            }
          }
        }
      }

      # Adaptive 2 starts with NB model, if it fails, it switches ZINB
      if (method == 'Adaptive2') {
        error1 <- try(obj1 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...))

        if (class(error1) == 'try-error' | obj1$pv == 0) {
          error2 <- try(obj2 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...))
          if (class(error2) != 'try-error') {
            obj <- obj2
            error <- error2
          } else {
            error <- error2
          }
        } else {
          error <- error1
          obj <- obj1
        }

      }

      # Random Effects model
      # ZINB, B, NB, Adpative1 is implemented based on glmmADMB

      # Set P value NA for those not makes sense
      if (class(error) == "try-error" | abs(obj$lfc) > 100) {
        obj$pv <- obj$lfc <- obj$lfc.lci <- obj$lfc.uci <- obj$method <- NA
      }


      pv.vec <- rbind(pv.vec, obj$pv)
      fc.vec <- rbind(fc.vec, obj$lfc)
      fc.lc.vec <- rbind(fc.lc.vec, obj$lfc.lci)
      fc.uc.vec <- rbind(fc.uc.vec, obj$lfc.uci)
      met.vec <- rbind(met.vec, obj$method)

    }
    cat('\n')

    qv.vec <- matrix(p.adjust(pv.vec[, 1], 'fdr'), ncol=1)

    rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(fc.uc.vec) <- rownames(fc.lc.vec) <- rownames(met.vec) <- rownames(prop)
    colnames(pv.vec) <- 'Pvalue'
    colnames(qv.vec) <- 'Qvalue'
    colnames(met.vec) <- 'Method'


    pv.list[[LOI]] <- pv.vec
    qv.list[[LOI]] <- qv.vec
    fc.list[[LOI]] <- fc.vec
    fc.lc.list[[LOI]] <- fc.lc.vec
    fc.uc.list[[LOI]] <- fc.uc.vec
    met.list[[LOI]] <- met.vec

    res <- cbind(pv.vec, qv.vec, fc.vec, fc.lc.vec, fc.uc.vec, met.vec)
    rownames(res) <- rownames(prop)
    write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))

    if (mt.method == 'fdr') {
      res.final <- rbind(res.final, res[as.numeric(res[, 'Qvalue']) <= cutoff, , drop=F])
    }
    if (mt.method == 'raw') {
      res.final <- rbind(res.final, res[ as.numeric(res[, 'Pvalue']) <= cutoff, , drop=F])
    }
  }

  if (!is.null(res.final)) {
    colnames(res.final) <- colnames(res)
    res.final <- res.final[rowSums(is.na(res.final)) == 0, , drop=F]
    write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
  }
  return(list(pv.list=pv.list, qv.list=qv.list, fc.list=fc.list, fc.uc.list=fc.uc.list, fc.lc.list=fc.lc.list, met.list=met.list))
}


# This function for nonparametric/permutaiton method Rev: 2017_02_16
# Add normalization method; Add transformation; Remove rarefaction
# (only output warnings); Rev: 2017_10_30 Support filtering based on
# coefficient of variation Rev: 2018_01_30 Add ct.min and handle NA
# size factor Rev: 2018_03_15 Add winsorization support for permutatin
# test - default is FALSE Rev: 2018_10_03 TSS size factor calculation
# is based on specific taxa level

#' Perform differential analysis
#'
#' @param data.obj Data object created by load_data()
#' @param grp.name Variable of interest
#' @param adj.name List of covariates / variables to adjust for
#' @param subject Variable indicating subject
#' @param taxa.levels Taxa levels to include in analysis
#' @param method Method to use for differential analysis (Default: "perm")
#' @param block.perm Perform block permutation (TRUE/FALSE, Default: FALSE)
#' @param perm.no Number of permutations to perform (Default: 999)
#' @param norm Normalization method (Default: "GMPR")
#' @param norm.level Normalization method (Default: "Species")
#' @param intersect.no (For GMPR) The minimum number of shared features between sample pair, where the ratio is calculated (Default: 4)
#' @param ct.min (For GMPR) The minimum number of counts required to calculate ratio (Default: 2)
#' @param transform Type of transformation to perform (Default: "sqrt")
#' @param winsor Perform Winsorization of data? (TRUE/FALSE, Default: FALSE)
#' @param winsor.qt Winsorization quartile (Default: 0.97)
#' @param prev Prevalence cutoff for OTUs (Default: 0.1)
#' @param minp Minimum proportion/abundance cutoff (Default: 0.002)
#' @param medianp Median proportion/abundance cutoff (Optional, can be used in place of \code{minp}. Default: NULL)
#' @param cv Coefficient of variance to use as OTU/ASV cutoff (Optional, can be used in place of \code{minp}. Default: NULL)
#' @param mt.method Method to use for multiple testing correction (Default: "fdr")
#' @param cutoff Q-value cutoff for significance testing after multiple testing correction
#' @param seed Random seed
#' @param ... Any additional parameters
#'
#' @return A dataframe containing the results of differential abundance testing
#' @export
#'
#' @examples
#' data("Constipation")
#' diff.obj <- perform_differential_analysis(data.obj, grp.name="Visit",
#'                                           taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
#'                                           method='perm', mt.method='fdr',
#'                                           block.perm=FALSE, norm = 'GMPR')
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
        prop0 <- prop0[matrixStats::rowMaxs(prop0) > minp, , drop = FALSE]
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
        prop0 <- prop0[matrixStats::rowMaxs(prop0) > minp, , drop = FALSE]
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


#' Visualize differential analysis results
#'
#' @param data.obj Data object created by load_data()
#' @param diff.obj Differential abundance object created by perform_differential_analysis()
#' @param grp.name Variable of interest
#' @param strata Variable indicating strata (Optional, default: NULL)
#' @param test Type of test, either "Para" or "Nonpara" (Default: "Nonpara")
#' @param mt.method Type of multiple testing method used when creating \code{diff.obj}
#' @param scale Scaling method used when creating \code{diff.obj}
#' @param cutoff Q-value cutoff used when creating \code{diff.obj}
#' @param taxa.levels Taxa levels to visualize (Default includes "Phylum", "Family", and "Genus")
#' @param ord Order by mean proportion, TRUE/FALSE (Default: TRUE)
#' @param eff.type Effect type (Default: "logP")
#' @param indivplot Generate plots for individual samples, TRUE/FALSE (Default: TRUE)
#' @param subject Variable indicating subject
#'
#' @return A list of differential abundance plots (boxplots, barplots, heatmaps, biplots, effect size plots)
#' @export
#'
#' @examples
#' data("Constipation")
#' diff.obj <- perform_differential_analysis(data.obj, grp.name="Visit",
#'                                           taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
#'                                           method='perm', mt.method='fdr',
#'                                           block.perm=FALSE, norm = 'GMPR')
#' diff_vis <- visualize_differential_analysis(data.obj, diff.obj, grp.name="Visit", taxa.levels=c('Phylum', 'Family', 'Genus'), mt.method='fdr', cutoff=0.2049)
visualize_differential_analysis <- function (data.obj, diff.obj,  grp.name=NULL, strata=NULL, test='Nonpara', mt.method='fdr', scale='sqrt', cutoff=0.15,
                                             taxa.levels=c('Phylum', 'Family', 'Genus'), ord=TRUE, eff.type='logP', indivplot=TRUE, subject=NULL) {

  # uniquefy names
  # For backward compatibility. Newer version will not need this and below. The old version has 'unclassified' which leads to duplicate names.
  # Newer version has 'Unclassified'. Case difference.

  # Check whether there is name duplication
  check.names <- NULL
  results <- list()
  obj0 <- diff.obj[[1]]
  for (level in names(obj0)) {
    obj <- obj0[[level]]
    # rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
    check.names <- c(check.names, rownames(obj))
  }

  if (sum(table(check.names) >= 2)) {
    data.obj <- uniquefy_taxa_names(data.obj)

    for (name1 in names(diff.obj)) {
      obj0 <- diff.obj[[name1]]
      for (level in names(obj0)) {
        obj <- obj0[[level]]
        # rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
        rownames(obj) <- paste0(rownames(obj), substr(level, 1, 1))
        obj0[[level]] <- obj
      }
      diff.obj[[name1]] <- obj0
    }

  }

  fc.list <- diff.obj$fc.list
  qv.list <- diff.obj$qv.list
  pv.list <- diff.obj$pv.list
  if (test == 'Para') {
    fc.lc.list <- diff.obj$fc.lc.list
    fc.uc.list <- diff.obj$fc.uc.list
  }
  df <- data.obj$meta.dat
  grp <- df[, grp.name]

  ind <- !is.na(grp)
  data.obj <- subset_data(data.obj, ind)
  grp <- grp[ind]
  df <- df[ind, ]

  prop <- NULL
  eff <- eff.lc <- eff.uc <- NULL
  taxa.names <- NULL
  if (is.null(taxa.levels)) {
    LOIs <- names(qv.list)
  } else {
    LOIs <- taxa.levels
    if (sum(!(taxa.levels %in% names(qv.list)))) {
      stop('Taxa levels are not contained in differential abundance analysis results!\n')
    }
  }
  for (LOI in LOIs) {
    pv.vec <- pv.list[[LOI]]
    fc.vec <- fc.list[[LOI]]
    #qv.vec <- qvalue(pv.vec[, 1])$qvalues
    qv.vec <- qv.list[[LOI]]

    if (test == 'Para') {
      fc.lc.vec <- fc.lc.list[[LOI]]
      fc.uc.vec <- fc.uc.list[[LOI]]
    }

    if (mt.method == 'fdr') {
      taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
      taxa.name <- taxa.name[!is.na(taxa.name)]
    }

    if (mt.method == 'raw') {
      taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
      taxa.name <- taxa.name[!is.na(taxa.name)]
    }

    if (length(taxa.name) != 0) {
      prop0 <- data.obj$abund.list[[LOI]]
      prop0 <- t(t(prop0) / colSums(prop0))
      prop0 <-  prop0[taxa.name, , drop=F]
      if (ord == TRUE) {
        prop0 <- prop0[rev(order(rowMeans(prop0))), , drop=F]
      }
      prop <- rbind(prop, prop0)
      # currently using fold change
      if (test == 'Para') {
        eff <- rbind(eff, fc.vec[taxa.name, , drop=F])
        eff.lc <- rbind(eff.lc, fc.lc.vec[taxa.name, , drop=F])
        eff.uc <- rbind(eff.uc, fc.uc.vec[taxa.name, , drop=F])
      } else {
        if (eff.type == 'LFC') {
          eff <- c(eff, fc.vec[taxa.name, ])
        }
        if (eff.type == 'logP') {
          eff <- c(eff, sign(fc.vec[taxa.name, ]) * (-log10(pv.vec[taxa.name, ])))
        }
      }
      taxa.names <- c(taxa.names, taxa.name)
    }
  }
  results$taxa.names <- taxa.names
  if (length(taxa.names) == 0) {
    cat('No differential taxa! \n')
  } else {
    if (length(taxa.names) >= 2) {

      if (!is.null(grp.name)) {
        results$barplot_aggregate <- taxa_barplot_aggregate(prop, df, grp.name, strata, scale)
        results$boxplot_aggregate <- taxa_boxplot_aggregate(prop, df, grp.name, strata, scale)
      }

      # currently fold change
      if (test == 'Para') {
        rownames(eff) <- rownames(eff.lc) <- rownames(eff.uc) <- taxa.names
        for (k in 1:ncol(eff)) {
          fold.dat.plot1 <- data.frame(Estimate=eff[, k], LCI=eff.lc[, k], UCI=eff.uc[, k], IV=taxa.names)
          results$effect_size <- plot_effect_size2(fold.dat.plot1)
        }
      } else {
        if (!is.na(eff[1])) {
          names(eff) <- taxa.names
          eff <- eff[!is.na(eff) & is.finite(eff)]
          eff <- sort(eff)
          taxa.names2 <- names(eff)

          if (eff.type == 'LFC') {
            results$effect_size <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Log2 fold change')
          }
          if (eff.type == 'Spearman') {
            results$effect_size <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Spearman correlation')
          }
          if (eff.type == 'logP') {
            results$effect_size <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='-log10(P)')
          }
        }
      }

      # create heatmp
      #			taxa.names2 <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]
      if (!is.null(grp.name)) {

        results$prop_heatmap <- generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=c(grp.name, strata),  ann=NULL)
        results$rank_heatmap <-generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=c(grp.name, strata), data.type='R', ann=NULL)
        try(
          results$biplot <- generate_taxa_biplot(data.obj, taxa=taxa.names, trans='sqrt', grp.name, ann=NULL, varname.size = 1.5)
        )
      }
    }
    if (!is.null(grp.name)) {
      # Individual plots
      if (indivplot == TRUE) {
        results$taxa_boxplot <- generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, ann=NULL)
        if (!is.null(subject)) {
          generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, subject=subject, ann=NULL)
        }
        results$taxa_boxplot_binary <- generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, scale='binary', ann=NULL)
        results$taxa_barplot <- generate_taxa_barplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, ann=NULL)
      }

    }
  }
  #if(length(intersect(taxa.levels, c('Phylum', 'Family', 'Genus')))){
  #  results$cladogram <- make_cladogram(data.obj=data.obj, diff.obj=diff.obj, grp.name=grp.name, mt.method=mt.method, cutoff=cutoff)
  #}
  return(results)
}

# Rev: 2016_11_25 Uniqufy
# Rev: 2016_04_18 strata

OLD.visualize_differential_analysis <- function (data.obj, diff.obj,  grp.name=NULL, strata=NULL, test='Nonpara', mt.method='fdr', scale='sqrt', cutoff=0.15,
                                             taxa.levels=c('Phylum', 'Family', 'Genus'), ord=TRUE, eff.type='logP', indivplot=TRUE, colFnsC=NULL, colFnsF=NULL, subject=NULL,
                                             xsize=10, ann='', hei1=NULL, wid1=NULL, hei2=NULL, wid2=NULL) {

  # uniquefy names
  # For backward compatibility. Newer version will not need this and below. The old version has 'unclassified' which leads to duplicate names.
  # Newer version has 'Unclassified'. Case difference.

  # Check whether there is name duplication
  check.names <- NULL

  obj0 <- diff.obj[[1]]
  for (level in names(obj0)) {
    obj <- obj0[[level]]
    # rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
    check.names <- c(check.names, rownames(obj))
  }

  if (sum(table(check.names) >= 2)) {
    data.obj <- uniquefy_taxa_names(data.obj)

    for (name1 in names(diff.obj)) {
      obj0 <- diff.obj[[name1]]
      for (level in names(obj0)) {
        obj <- obj0[[level]]
        # rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
        rownames(obj) <- paste0(rownames(obj), substr(level, 1, 1))
        obj0[[level]] <- obj
      }
      diff.obj[[name1]] <- obj0
    }

  }

  #pv.list <- diff.obj$pv.list
  fc.list <- diff.obj$fc.list

  if (mt.method == 'fdr.perm') {
    qv.list <- diff.obj$qv.perm.list
  } else {
    qv.list <- diff.obj$qv.list
  }

  pv.list <- diff.obj$pv.list
  if (test == 'Para') {
    fc.lc.list <- diff.obj$fc.lc.list
    fc.uc.list <- diff.obj$fc.uc.list
  }

  R2.list <- diff.obj$R2.list

  df <- data.obj$meta.dat
  grp <- df[, grp.name]

  ind <- !is.na(grp)
  data.obj <- subset_data(data.obj, ind)
  grp <- grp[ind]
  df <- df[ind, ]

  prop <- NULL
  eff <- eff.lc <- eff.uc <- R2 <- NULL
  taxa.names <- NULL
  if (is.null(taxa.levels)) {
    LOIs <- names(qv.list)
  } else {
    LOIs <- taxa.levels
    if (sum(!(taxa.levels %in% names(qv.list)))) {
      stop('Taxa levels are not contained in differential abundance analysis results!\n')
    }
  }
  for (LOI in LOIs) {
    pv.vec <- pv.list[[LOI]]
    fc.vec <- fc.list[[LOI]]
    #qv.vec <- qvalue(pv.vec[, 1])$qvalues
    qv.vec <- qv.list[[LOI]]

    if (!is.null(R2.list)) R2.vec <- R2.list[[LOI]]

    if (test == 'Para') {
      fc.lc.vec <- fc.lc.list[[LOI]]
      fc.uc.vec <- fc.uc.list[[LOI]]
    }

    if (mt.method %in% c('fdr', 'fdr.perm')) {
      taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
      taxa.name <- taxa.name[!is.na(taxa.name)]
    }

    if (mt.method == 'raw') {
      taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
      taxa.name <- taxa.name[!is.na(taxa.name)]
    }

    if (length(taxa.name) != 0) {
      prop0 <- data.obj$abund.list[[LOI]]
      prop0 <- t(t(prop0) / colSums(prop0))
      prop0 <-  prop0[taxa.name, , drop=F]
      if (ord == TRUE) {
        prop0 <- prop0[rev(order(rowMeans(prop0))), , drop=F]
      }
      prop <- rbind(prop, prop0)
      # currently using fold change
      if (test == 'Para') {
        eff <- rbind(eff, fc.vec[taxa.name, , drop=F])
        eff.lc <- rbind(eff.lc, fc.lc.vec[taxa.name, , drop=F])
        eff.uc <- rbind(eff.uc, fc.uc.vec[taxa.name, , drop=F])
      } else {
        if (eff.type %in% c('LFC', 'Spearman')) {
          eff <- c(eff, fc.vec[taxa.name, ])
        }
        if (eff.type == 'logP') {
          eff <- c(eff, sign(fc.vec[taxa.name, ]) * (-log10(pv.vec[taxa.name, ])))
        }
        if (!is.null(R2.list)) R2 <- c(R2, R2.vec[taxa.name, ])

      }
      taxa.names <- c(taxa.names, taxa.name)
    }
  }

  if (length(taxa.names) == 0) {
    cat('No differential taxa! \n')
  } else {
    if (length(taxa.names) >= 2) {
      if (is.null(wid1) | is.null(hei1)) {
        wid1 <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
        hei1 <- 7
      }

      if (!is.null(grp.name)) {
        obj <- taxa_barplot_aggregate(prop, df, grp.name, strata, scale, xsize)
        pdf(paste("Taxa_DifferentialAbundance_AbundanceBarplot", scale, mt.method, cutoff, ann, ".pdf", sep="_"), height=hei1, width=wid1)
        print(obj)
        dev.off()

        obj1 <- taxa_boxplot_aggregate(prop, df, grp.name, strata, scale, xsize)
        pdf(paste("Taxa_DifferentialAbundance_AbundanceBoxplot", scale, mt.method, cutoff, ann, ".pdf", sep="_"), height=hei1, width=wid1)
        print(obj1)
        dev.off()

      }

      # currently fold change
      if (test == 'Para') {
        rownames(eff) <- rownames(eff.lc) <- rownames(eff.uc) <- taxa.names
        for (k in 1:ncol(eff)) {
          fold.dat.plot1 <- data.frame(Estimate=eff[, k], LCI=eff.lc[, k], UCI=eff.uc[, k], IV=taxa.names)
          obj <- plot_effect_size2(fold.dat.plot1)
          pdf(paste("Taxa_DifferentialAbundance_EffectBarplot", colnames(eff)[k], mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
          print(obj)
          dev.off()
        }
      } else {
        if (!is.na(eff[1])) {
          names(eff) <- taxa.names
          eff <- eff[!is.na(eff) & is.finite(eff)]
          eff <- sort(eff)
          taxa.names2 <- names(eff)
          if (is.null(wid2) | is.null(hei2)) {
            hei2 <- 4 + length(taxa.names2) / 20 * 3
            wid2 <- 6
          }
          if (eff.type == 'LFC') {
            obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Log2 fold change')
            pdf(paste("Taxa_DifferentialAbundance_LFCBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
            print(obj)
            dev.off()
          }

          if (eff.type == 'Spearman') {
            obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Spearman correlation')
            pdf(paste("Taxa_DifferentialAbundance_SpearmanBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
            print(obj)
            dev.off()
          }

          if (eff.type == 'logP') {
            obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='-log10(P)')
            pdf(paste("Taxa_DifferentialAbundance_logPBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
            print(obj)
            dev.off()
          }

          if (!is.null(R2.list))   {
            R2 <- sort(R2, decr=TRUE)
            taxa <- names(R2)
            taxa <- factor(taxa, levels=taxa)
            df <- data.frame(taxa = taxa, R2 = R2)
            obj <- ggplot(df, aes(taxa, R2)) + geom_col(fill='steelblue', col='gray', size=0.25) + ylab(expression(R^2)) + xlab('') +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
            pdf(paste("Taxa_DifferentialAbundance_R2_Barplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=wid2, width=hei2)
            print(obj)
            dev.off()

          }


        }
      }

      # create heatmp
      #			taxa.names2 <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]
      if (!is.null(grp.name)) {

        generate_taxa_heatmap(data.obj, taxa.levels='All', sam.ord=order(grp), taxa=taxa.names, meta.info=c(grp.name, strata), Colv=F, dendrogram='row',
                              ann=paste0(mt.method, '_', cutoff, '_', ann, '_Unclustered'), colFnsC=colFnsC, colFnsF=colFnsF)
        generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=c(grp.name, strata),  ann=paste0(mt.method, '_', cutoff, '_', ann), colFnsC=colFnsC, colFnsF=colFnsF)
        generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=c(grp.name, strata), data.type='R', ann=paste0(mt.method, '_', cutoff, '_Rank_', ann), colFnsC=colFnsC, colFnsF=colFnsF)
        try(
          generate_taxa_biplot(data.obj, taxa=taxa.names, trans='sqrt', grp.name, ann=paste0(mt.method, '_', cutoff, '_', ann), varname.size = 1.5)
        )
      }
    }
    if (!is.null(grp.name)) {
      # Individual plots
      if (indivplot == TRUE) {
        generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, ann=paste0(mt.method, '_', cutoff, '_', ann))
        if (!is.null(subject)) {
          generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, subject=subject, ann=paste0(mt.method, '_', cutoff, '_', ann, '_Paired'))
        }
        generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, scale='binary', ann=paste0(mt.method, '_', cutoff, '_', ann))
        generate_taxa_barplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, ann=paste0(mt.method, '_', cutoff, '_', ann))
      }

    }
  }
}


# Rev: 2016_12_07 when the names of some differential taxa are not consistent
# Rev: 2017_02_21  Color by log fold change m.list ->> fc.list
# To-do-list:  size by effects, accomodate OTU level
create_lefse_format <- function(data.obj, diff.obj, grp.name, mt.method='fdr', cutoff=0.15, prev=0.1, minp=0.002,
                                lefse.dir="/data2/microbiome/jeff/tools/nsegata-lefse/", ann="") {

  # To be improved - currently no effect size
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  levels(grp) <- paste0(1:nlevels(grp), levels(grp))

  qv.list <- diff.obj$qv.list
  pv.list <- diff.obj$pv.list
  fc.list <- diff.obj$fc.list

  otu.name.12 <- data.obj$otu.name
  otu.tab.12 <- data.obj$otu.tab
  tax.family.a <- NULL
  for (i in 1:6) {
    tax.family <- apply(otu.name.12[, 1:i, drop=F], 1, paste, collapse=".")
    phlan.tab <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
    rownames(phlan.tab) <- phlan.tab[, 1]
    phlan.tab <- as.matrix(phlan.tab [, -1, drop=F])
    phlan.tab <- t(t(phlan.tab) / colSums(phlan.tab))
    phlan.tab <- phlan.tab[rowSums(phlan.tab!=0) > prev*ncol(phlan.tab) & rowMaxs(phlan.tab) > minp, , drop=F]
    tax.family.a <- c(tax.family.a, rownames(phlan.tab))
    #tax.family.a <- c(tax.family.a, unique(tax.family))
  }
  alias.a <- sapply(strsplit(tax.family.a, "\\."), function(x) {
    if (length(x) >= 2) {
      if (length(x) == 2) {
        return(x[2])
      } else {
        return(paste0(x[2], ";", x[length(x)]))
      }
    } else {
      return(x)
    }
  })

  taxa.names <- NULL
  abundant.grp.names <- NULL

  for (LOI in setdiff(names(qv.list), 'Species')) {
    qv.vec <- qv.list[[LOI]]
    pv.vec <- pv.list[[LOI]]
    #qv.vec <- qvalue(pv.vec[, 1])$qvalues
    fc.vec <- fc.list[[LOI]]
    if (mt.method == 'fdr') {
      taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
      #	abundant.grp.name <- apply(m.vec, 1, function(x) levels(grp)[which.max(x)])[qv.vec <= cutoff]
      abundant.grp.name <- apply(fc.vec, 1, function(x) levels(grp)[as.numeric(x >= 0) + 1])[qv.vec <= cutoff]
    }

    if (mt.method == 'raw') {
      taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
      #	abundant.grp.name <- apply(m.vec, 1, function(x) levels(grp)[which.max(x)])[pv.vec <= cutoff]
      abundant.grp.name <- apply(fc.vec, 1, function(x) levels(grp)[as.numeric(x >= 0) + 1])[pv.vec <= cutoff]
    }

    if (length(taxa.name) != 0) {
      taxa.names <- c(taxa.names, taxa.name)
      abundant.grp.names <- c(abundant.grp.names, abundant.grp.name)
    }
  }

  # remove 'unclassified'
  abundant.grp.names <- abundant.grp.names[!grepl('unclassified', taxa.names, ignore.case=T)]
  taxa.names <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]

  # Rev: 2016_12_07
  #	taxa.names <- intersect(taxa.names, alias.a)

  ind <- match(taxa.names, alias.a)
  na.ind <- !is.na(ind)

  phlan <- cbind(tax.family.a, '1.5', '\t', '\t', '-')

  phlan[ind[na.ind], 2:5] <- cbind('3.5',  abundant.grp.names[na.ind], '3.0', '-')
  write.table(phlan,  paste0('Lefse.LDA.', ann, '.txt'), row.names=F, col.names=F, quote=F, sep='\t')
  #	cmd1 <- paste0('python ', lefse.dir,  'plot_cladogram.py Lefse.LDA.txt Taxa_DifferentialAbundance_', mt.method, cutoff, '_Cladogram_', ann, '.pdf', ' --format pdf')
  #	cmd2 <- paste0("python ", lefse.dir, 'plot_res.py Lefse.LDA.txt Taxa_DifferentialAbundance_', mt.method, cutoff, '_LDA_', ann, '.pdf --format pdf')
  #	system(cmd1)
  #	system(cmd2)
}

perform_lefse_analysis <- function (data.obj,  grp.name, sub.grp.name=NULL, prev=0.1, minp=0.002,
                                    lefse.dir="/data2/microbiome/jeff/tools/nsegata-lefse/", ann="") {
  otu.name.12 <- data.obj$otu.name
  otu.tab.12 <- data.obj$otu.tab
  meta.dat <- data.obj$meta.dat

  levels(meta.dat[, grp.name]) <- paste0(1:nlevels(meta.dat[, grp.name]), levels(meta.dat[, grp.name]))

  tax.family <- apply(otu.name.12, 1, paste, collapse="|")
  lefse.tab <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
  rownames(lefse.tab) <- lefse.tab[, 1]
  lefse.tab <- as.matrix(lefse.tab [, -1])

  lefse.tab <- t(t(lefse.tab) / colSums(lefse.tab))

  lefse.tab <- lefse.tab[rowMaxs(lefse.tab) > minp & rowSums(lefse.tab!=0) > prev*ncol(lefse.tab), , drop=FALSE]

  if (is.null(sub.grp.name)) {
    header <- rbind(class=as.character(meta.dat[colnames(lefse.tab), grp.name]),
                    id=colnames(lefse.tab))
  } else {
    header <- rbind(class=as.character(meta.dat[colnames(lefse.tab), grp.name]),
                    subclass=as.character(meta.dat[colnames(lefse.tab), sub.grp.name]),
                    id=colnames(lefse.tab))
  }


  lefse.tab <- rbind(header, lefse.tab)
  write.table(lefse.tab, "lefse.txt", sep="\t", col.names=F, quote=F)

  #	system(paste0('mkdir LefSe_', ann))
  #	output <- paste0('LefSe_', ann, '/')
  #	if (is.null(sub.grp.name)) {
  #		cmd1 <- paste0("python ", lefse.dir, "format_input.py lefse.txt ", output, "temp.in ", "-c 1 -u 2 -o 1000000")
  #	} else {
  #		cmd1 <- paste0("python ", lefse.dir, "format_input.py lefse.txt ", output, "temp.in ", "-c 1 -s 2 -u 3 -o 1000000")
  #	}
  #
  #	cmd2 <- paste0("python ", lefse.dir, "run_lefse.py ", output, "temp.in ",  output, "lda.res")
  #	cmd3 <- paste0("python ", lefse.dir, "plot_res.py ",  output, "lda.res ", output, "lda.pdf --format pdf")
  #	cmd4 <- paste0("python ", lefse.dir, "plot_cladogram.py ",  output, "lda.res ", output, "cladogram.pdf --format pdf --labeled_stop_lev 6 --abrv_stop_lev  6")
  #
  #	system(cmd1)
  #	system(cmd2)
  #	system(cmd3)
  #	system(cmd4)
}


plot.Boruta2 <- function (x, colCode = c("green", "yellow", "red", "blue"), sort = TRUE,
                          whichShadow = c(TRUE, TRUE, TRUE), col = NULL, xlab = "Attributes", ids=NULL,
                          ylab = "Importance", ...)
{
  if (class(x) != "Boruta")
    stop("This function needs Boruta object as an argument.")
  lz <- lapply(1:ncol(x$ImpHistory), function(i) x$ImpHistory[is.finite(x$ImpHistory[, i]), i])
  names(lz) <- colnames(x$ImpHistory)
  numShadow <- sum(whichShadow)
  lz <- lz[c(rep(TRUE, length(x$finalDecision)), whichShadow)]
  col <- Boruta:::generateCol(x, colCode, col, numShadow)
  if (sort) {
    ii <- order(sapply(lz, median))
    lz <- lz[ii]
    col <- col[ii]
  }
  names(lz) <- gsub("^X",  "", names(lz))
  if (is.null(ids)) {
    len <- sum(x$finalDecision %in% c('Confirmed', 'Tentative'))
    ind <- (length(lz) - len) : length(lz)
  } else {
    ind <- match(ids, names(lz))
  }

  boxplot(lz[ind], xlab = xlab, ylab = ylab, col = col[ind], ...)
  invisible(x)
  names(lz[ind])
}



createROC <- function (pv.list, lab.list, pos.lab='1', file.name='ROC.pdf', width = 6, height = 6) {
  require(ROCR)
  n <- length(pv.list)
  aucs <- numeric(n)
  names(aucs) <- names(pv.list)

  #	cols <- scales::hue_pal()(n)
  cols <- rep(c('red', 'blue', 'orange', 'cyan', 'purple'), ceiling(n/5))[1:n]
  ltys <- rep(c(1, 2), ceiling(n/2))[1:n]
  pdf(file.name, height=6, width=6)
  for (i in 1:n) {

    cat("*")
    pv.mat <- pv.list[[i]]
    lab.mat <- lab.list[[i]]

    pred <- prediction(pv.mat, lab.mat==pos.lab)
    perf <- performance(pred, "tpr", "fpr")
    aucs[i] <- mean(unlist(performance(pred, 'auc')@y.values))
    plot(perf, avg="threshold", spread.estimate = 'stddev', spread.scale = 2,  col=cols[i], lty=ltys[i], lwd=1.5,  add=ifelse(i==1, FALSE, TRUE),  main='ROC curve')

  }
  abline(0, 1, col='black')
  legend("right", legend=paste0(names(pv.list), "(AUC:", round(aucs, 3), ")"), col=cols, lty=ltys, lwd=2,  bty="n")
  dev.off()
}

#invlogit <- function(x) {
#	rbinom(length(x), 1, exp(x) / (1 + exp(x)))
#}
#pv.list <- list(x1=rnorm(100), x2=rnorm(100))
#lab.list <- list(x1=invlogit(pv.list[['x1']]), x2=invlogit(pv.list[['x2']]))
#createROC(pv.list, lab.list)

# Rev: 2017_02_17 supply aug.var to decision tree
# Rev: 2017_11_08 Change formula
predictionRF <- function (data.obj,  resp.name, formula=NULL, taxa.level='Species', binary=FALSE, prev=0.1, minp=0.002, B=50, seed=123,
                          boruta.level=c('Confirmed', 'Tentative'), ann='',...) {

  #sink(paste0("Taxa_RandomForest_", taxa.level, ".txt"))
  date()
  response <- data.obj$meta.dat[, resp.name]

  if (!is.null(formula)) {
    # adj.var <- unlist(strsplit(unlist(strsplit(formula, '\\s*~\\s*'))[2], '\\s*\\+\\s*'))
    # adj <- data.obj$meta.dat[, adj.var, drop=F]
    adj <- model.matrix(as.formula(formula), data.obj$meta.dat)
    adj.var <- colnames(adj)
  }

  if (taxa.level == 'Species') {
    if (taxa.level %in% names(data.obj$abund.list)) {
      ct <- data.obj$abund.list[[taxa.level]]
    } else {
      # Accomodate different version
      ct <- data.obj$otu.tab
      rownames(ct) <- paste0("OTU", rownames(ct), ":", data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
      data.obj$abund.list[['Species']] <- ct
    }
  } else {
    ct <- data.obj$abund.list[[taxa.level]]
  }

  prop <- t(t(ct) / colSums(ct))
  prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
  prop <- t(prop)

  if (binary == TRUE) {
    prop <- (prop != 0)
  }

  original.names <- colnames(prop)
  set.seed(seed)

  if (is.null(formula)) {
    performance <- matrix(0, nrow=B, ncol=2)
    colnames(performance) <- c("RF_M","Guess")
    roc.list <- list(RF_M=NULL)
    lab.list <- roc.list
  } else {
    performance <- matrix(0, nrow=B, ncol=4)
    colnames(performance) <- c("RF_M", "RF_CF", "RF_M+CF", "Guess")
    roc.list <- list("RF_M"=NULL, "RF_CF"=NULL, "RF_M+CF"=NULL)
    lab.list <- roc.list
  }


  colnames(prop) <- gsub(";", "_", colnames(prop))
  colnames(prop) <- gsub(":", "_", colnames(prop))
  colnames(prop) <- gsub("-", "_", colnames(prop))
  colnames(prop) <- gsub("\\.", "_", colnames(prop))
  names(original.names) <- colnames(prop)
  if (!is.null(formula)) {
    names(adj.var) <- adj.var
    original.names <- c(original.names, adj.var)
  }


  if (!is.null(formula)) {
    padj <- cbind(prop, adj)
  } else {
    padj <- prop
  }

  I <- nrow(prop)
  cat('Begin to bootstrap ...\n')
  for(b in 1:B){
    if (b %% 50 == 0) cat (".")
    err <- try({
      bsample <- sample(1:I, I, replace=T)
      if (is.factor(response)) {
        rf1 <- randomForest(x=prop[bsample, ], y=response[bsample], xtest=prop[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
        performance[b, "RF_M"] <- mean(response[-bsample]!= rf1$test$predicted)
        performance[b, "Guess"]<- mean(response[-bsample] != levels(response)[which.max(tabulate(response[bsample]))])
        roc.list[['RF_M']] <- c(roc.list[['RF_M']], rf1$test$vote[, 1])
        lab.list[['RF_M']] <- c(lab.list[['RF_M']], response[-bsample])
        if (!is.null(formula)) {
          rf2 <- randomForest(x=adj[bsample, ], y=response[bsample], xtest=adj[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
          rf3 <- randomForest(x=padj[bsample, ], y=response[bsample], xtest=padj[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
          performance[b, "RF_CF"] <- mean(response[-bsample]!= rf2$test$predicted)
          performance[b, "RF_M+CF"] <- mean(response[-bsample]!= rf3$test$predicted)
          roc.list[['RF_CF']] <- c(roc.list[['RF_CF']], rf2$test$vote[, 1])
          roc.list[['RF_M+CF']] <- c(roc.list[['RF_M+CF']], rf3$test$vote[, 1])
          lab.list[['RF_CF']] <- c(lab.list[['RF_CF']], response[-bsample])
          lab.list[['RF_M+CF']] <- c(lab.list[['RF_M+CF']], response[-bsample])
        }

      } else {
        rf1 <- randomForest(x=prop[bsample, ], y=response[bsample], xtest=prop[-bsample, ], ytest=response[-bsample], ...)
        performance[b, "RF_M"] <- mean((response[-bsample] - rf1$test$predicted)^2)
        performance[b, "Guess"]<- mean((response[-bsample] - mean(response[bsample]))^2)

        if (!is.null(formula)) {

          rf2 <- randomForest(x=adj[bsample, ], y=response[bsample], xtest=adj[-bsample, ], ytest=response[-bsample],  ...)
          rf3 <- randomForest(x=padj[bsample, ], y=response[bsample], xtest=padj[-bsample, ], ytest=response[-bsample],  ...)
          performance[b, "RF_CF"] <- mean((response[-bsample] - rf2$test$predicted)^2)
          performance[b, "RF_M+CF"] <- mean((response[-bsample]- rf3$test$predicted)^2)
        }
      }
    })
    if (inherits(err, "try-error")) {
      next
    }

  }
  cat("\n")

  sink(paste0("Taxa_Random_forest_FriedmanTest_", taxa.level, '_', ann, ".txt"))
  if (is.null(formula)) {
    cat("Fridman.test p value: ", friedman.test(performance)$p.value, "\n")
  } else {
    cat("Fridman.test p value (M+CF vs CF): ", friedman.test(performance[, c('RF_M+CF', 'RF_CF')])$p.value, "\n")
  }
  sink()

  pdf(paste0("Taxa_Random_forest_misclassification_barplot_", taxa.level, '_', ann, ".pdf"), height=5, width=4)
  if (!is.null(formula)) {
    performance2 <- performance[, c('Guess', 'RF_CF', 'RF_M', 'RF_M+CF')]
  } else {
    performance2 <- performance
  }
  if (is.factor(response)) {
    boxplot(performance2, col="#4DAF4A", ylab='Classification error', las=2)
  } else {
    boxplot(performance2, col="#4DAF4A", ylab='PMSE', las=2)
  }
  dev.off()

  # ROC curve - only compare to the first level of the factor
  if (is.factor(response)) {
    lab.list <- lapply(lab.list, function(x) {
      y <- factor(x)
      levels(y) <- levels(response)
      y
    }
    )
    createROC(roc.list, lab.list, pos.lab=levels(response)[1],
              file.name=paste0("Taxa_Random_forest_ROC_", taxa.level, '_', ann, ".pdf"))
  }

  if (is.factor(response)) {
    rf <- randomForest(x=padj, y=response, sampsize=table(response), importance=T, ...)
  } else {
    rf <- randomForest(x=padj, y=response,  importance=T, ...)
  }

  importance <- as.data.frame(rf$importance)

  if (is.factor(response)) {
    write.csv(importance[rev(order(importance$MeanDecreaseAccuracy)), ],
              paste0("Taxa_RandomForest_Ranking_MeanDecreaseAccuracy_", taxa.level,'_', ann, ".csv"))
    write.csv(importance[rev(order(importance$MeanDecreaseGini)), ],
              paste0("Taxa_RandomForest_Ranking_MeanDecreaseGini_", taxa.level, '_', ann, ".csv"))
  } else {
    write.csv(importance[rev(order(importance[, '%IncMSE'])), ],
              paste0("Taxa_RandomForest_Ranking_IncMSE_", taxa.level, '_', ann, ".csv"))
    write.csv(importance[rev(order(importance$IncNodePurity)), ],
              paste0("Taxa_RandomForest_Ranking_IncNodePurity_", taxa.level, '_', ann, ".csv"))
  }


  # Boruta Feature Selection - All subset selection (can't solve the confounding problem)
  obj.Boruta <- Boruta(padj, response, doTrace = 2)
  write.csv(obj.Boruta$finalDecision, paste0("Taxa_Random_forest_Boruta_Feature_Selection_", taxa.level, '_', ann, ".csv"))

  pdf(paste0("Taxa_Random_forest_Boruta_Feature_Selection_", taxa.level, '_', ann,  ".pdf"), height=6, width=10)
  par(mar=par('mar') + c(3, 0, 0, 0))
  plot(obj.Boruta, main = "Feature selection by Boruta", ylab="Importance z-score", lwd = 0.5, las = 3, xlab = "",
       cex=1 / (ncol(prop)/50), cex.axis=0.25*200/ncol(prop), yaxt='n')
  axis(2, cex.axis=1)
  dev.off()
  #sink()

  pdf(paste0("Taxa_Random_forest_Boruta_Feature_Selection_Significant_Only_", taxa.level, '_', ann, ".pdf"), height=6, width=6)
  par(mar=par('mar') + c(3, 0, 0, 0))
  otu.ids <- plot.Boruta2(obj.Boruta, main = "Feature selection by Boruta", lwd = 0.5, las = 3,
                          ylab="Importance z-score", xlab = "", cex.axis=1, yaxt='n')
  axis(2, cex.axis=1)
  dev.off()

  if ('Confirmed' %in% boruta.level) {
    taxa.names <- original.names[names(obj.Boruta$finalDecision)[obj.Boruta$finalDecision %in% c('Confirmed')]]
    if (!is.null(formula)) {
      taxa.names <- setdiff(taxa.names, adj.var)
      aug.var <- adj.var
    } else {
      aug.var <- NULL
    }

    if (length(taxa.names) > 0) {
      if (length(taxa.names) > 1) {

        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12,
                              margins=c(5, 15), ann=paste0('BorutaFeatures_P_Confirmed_', ann))
        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='R',
                              margins=c(5, 15), ann=paste0('BorutaFeatures_R_Comfirmed_', ann))
        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='B',
                              margins=c(5, 15), ann=paste0('BorutaFeatures_B_Comfirmed_', ann))
        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, Colv=F, dendrogram='row', sam.ord=order(response),
                              margins=c(5, 15), ann=paste0('BorutaFeatures_P_Confirmed_Unclusterded_', ann))
        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='R', Colv=F, dendrogram='row', sam.ord=order(response),
                              margins=c(5, 15), ann=paste0('BorutaFeatures_R_Comfirmed_Unclusterded_', ann))
        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='B', Colv=F, dendrogram='row', sam.ord=order(response),
                              margins=c(5, 15), ann=paste0('BorutaFeatures_B_Comfirmed_Unclusterded_', ann))
      }

      if (is.factor(response)) {
        generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Confirmed_', ann))
        generate_taxa_barplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        if (length(taxa.names) > 1) {
          generate_taxa_barplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
          generate_taxa_boxplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
          build.decision.tree(data.obj,  resp.name=resp.name, aug.var=aug.var, taxa.level=taxa.level, binary=binary, taxa=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        }
        gen <- data.obj$abund.list[[taxa.level]]
        gen <- t(t(gen) / colSums(gen))
        dat <- t(gen[taxa.names, ,drop=FALSE])
        colnames(dat) <- paste0('V', 1:ncol(dat))
        response2 <- factor(-as.numeric(response)+2)

        mylr <- function(formula, train, test){
          model <- randomForest(formula, data=train)
          x <- predict(model, newdata=test, type='prob')[, 2]
        }

        ACC <- Daim(response2 ~., model=mylr, data=as.data.frame(dat), labpos="1",
                    control=Daim.control(method="boot", number=100), cutoff="0.632+")
        pdf(paste0('BorutaFeatures_Confirmed_ROC_', taxa.level, '_0.632+', ann, '.pdf'), height=6, width=6)
        plot(ACC, main='Boruta taxa', method='0.632+', lwd=1.5)
        abline(0, 1, col='black')
        x <- ACC
        legend("bottomright", legend = paste("AUC:",
                                             formatC(Daim::auc(x$roc$sens632p, x$roc$spec632p),
                                                     digits = max(3, getOption("digits") - 3))),
               inset = 0.01)
        dev.off()
      } else {
        data.obj$meta.dat[, paste0(resp.name, '_b')] <- factor(data.obj$meta.dat[, resp.name] < median(data.obj$meta.dat[, resp.name]))
        levels(data.obj$meta.dat[, paste0(resp.name, '_b')]) <- c('High', 'Low')
        generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Confirmed_', ann))
        generate_taxa_barplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        if (length(taxa.names) > 1) {
          generate_taxa_barplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
          generate_taxa_boxplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        }
      }
    }
  }

  if ('Tentative' %in% boruta.level) {
    taxa.names <- original.names[names(obj.Boruta$finalDecision)[obj.Boruta$finalDecision %in% c('Confirmed',  'Tentative')]]
    if (!is.null(formula)) {
      taxa.names <- setdiff(taxa.names, adj.var)
      aug.var <- adj.var
    } else {
      aug.var <- NULL
    }

    if (length(taxa.names) > 0) {

      if (length(taxa.names) > 1)  {
        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12,
                              margins=c(5, 15), ann=paste0('BorutaFeatures_P_Tentative_', ann))
        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='R',
                              margins=c(5, 15), ann=paste0('BorutaFeatures_R_Tentative_', ann))
        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='B',
                              margins=c(5, 15), ann=paste0('BorutaFeatures_B_Tentative_', ann))

        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, Colv=F, dendrogram='row', sam.ord=order(response),
                              margins=c(5, 15), ann=paste0('BorutaFeatures_P_Tentative_Unclustered_', ann))
        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='R', Colv=F, dendrogram='row', sam.ord=order(response),
                              margins=c(5, 15), ann=paste0('BorutaFeatures_R_Tentative_Unclustered_', ann))
        generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='B', Colv=F, dendrogram='row', sam.ord=order(response),
                              margins=c(5, 15), ann=paste0('BorutaFeatures_B_Tentative_Unclustered_', ann))
      }

      if (is.factor(response)) {
        generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Tentative_', ann))
        generate_taxa_barplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        if (length(taxa.names) > 1) {
          generate_taxa_barplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
          generate_taxa_boxplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
          build.decision.tree(data.obj,  resp.name=resp.name, aug.var=aug.var, taxa.level=taxa.level, binary=binary, taxa=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        }

        gen <- data.obj$abund.list[[taxa.level]]
        gen <- t(t(gen) / colSums(gen))
        dat <- t(gen[taxa.names, ,drop=FALSE])
        colnames(dat) <- paste0('V', 1:ncol(dat))
        response2 <- factor(-as.numeric(response)+2)

        mylr <- function(formula, train, test){
          model <- randomForest(formula, data=train)
          x <- predict(model, newdata=test, type='prob')[, 2]
        }

        ACC <- Daim(response2 ~., model=mylr, data=as.data.frame(dat), labpos="1",
                    control=Daim.control(method="boot", number=100), cutoff="0.632+")
        pdf(paste0('BorutaFeatures_Tentative_ROC_', taxa.level, '_0.632+', ann, '.pdf'), height=6, width=6)
        plot(ACC, main='Boruta taxa', method='0.632+', lwd=2)
        abline(0, 1, col='black')
        x <- ACC
        legend("bottomright", legend = paste("AUC:",
                                             formatC(Daim::auc(x$roc$sens632p, x$roc$spec632p),
                                                     digits = max(3, getOption("digits") - 3))),
               inset = 0.01)
        dev.off()

      } else {
        data.obj$meta.dat[, paste0(resp.name, '_b')] <- factor(data.obj$meta.dat[, resp.name] < median(data.obj$meta.dat[, resp.name]))
        levels(data.obj$meta.dat[, paste0(resp.name, '_b')]) <- c('High', 'Low')
        generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Tentative_', ann))
        generate_taxa_barplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        if (length(taxa.names) > 1) {
          generate_taxa_barplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
          generate_taxa_boxplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        }
      }
    }
  }
}
