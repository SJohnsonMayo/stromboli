

# Rev: 2017_02_13 Individual label
# Rev: 2017_05_06 Add p value
generate_ordination <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                 grp.name, adj.name=NULL, emp.lev=NULL, indiv.lab=FALSE, indiv.lab.cex=0.5, is.lightcolor=TRUE,
                                 strata=NULL,  pc.separate, pca.method='cmd', ann=NULL, sub=NULL,
                                 clab=1.0, cex.pt=1.25, ellipse=T, cstar= 1, wid=5, hei=5, pdf=TRUE, is.pvalue=FALSE, ...) {
  # Implment strata
  # To be completed, add continuous case
  strata0 <- strata

  if (pdf) {
    if (is.null(ann)) {
      pdf(paste0('Beta_diversity_ordination_', pca.method, '_', grp.name, '.pdf'), width=wid, height=hei)
    } else {
      pdf(paste0('Beta_diversity_ordination_', pca.method, '_', ann, '.pdf'), width=wid, height=hei)
    }
  }

  df <- data.obj$meta.dat
  grp <- factor(df[, grp.name])

  if (!is.null(emp.lev)) {
    grp <- factor(grp, levels=c(setdiff(levels(grp), emp.lev), emp.lev))
  }

  if (!is.null(strata)) {
    strata <- factor(df[, strata])
  } else {
    strata <- factor(grp)
  }

  # To be revised
  darkcols <- hue_pal(l=40)(nlevels(grp))
  if (is.lightcolor) {
    lightcols <- hue_pal(c=45, l=80)(nlevels(grp))
  } else {
    lightcols <- darkcols
  }

  pchs <- rep(c(21, 22, 23, 24, 25), ceiling(nlevels(strata) / 5))[1:nlevels(strata)]

  for (dist.name in dist.names) {
    dist.temp <- dist.obj[[dist.name]]
    if (!is.null(adj.name)) {
      adj <- as.data.frame(df[, adj.name])
      obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
      dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
      dist.temp <- dist(dat2)
    }
    if (pca.method == 'cmd') {
      obj <- cmdscale(as.dist(dist.temp), k=2, eig=T)
      pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
      y <- cbind(obj$points[, 1], obj$points[, 2])

      xlab <- paste0('PC1(', pve[1], '%)')
      ylab <- paste0('PC2(', pve[2], '%)')
    }

    if (pca.method == 'nmds') {
      obj <- metaMDS(as.dist(dist.temp), k=2)
      y <- cbind(obj$points[, 1], obj$points[, 2])
      xlab <- 'NMDS1'
      ylab <- 'NMDS2'
    }

    if (pca.method == 'pls') {
      require(mixOmics)
      # Test
      obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
      obj <- plsda(obj, grp, ncomp=2)
      y <- cbind(obj$variates$X[, 1], obj$variates$X[, 2])
      xlab <- 'PLS1'
      ylab <- 'PLS2'
    }

    if (is.pvalue == TRUE & nlevels(factor(grp)) > 1) {
      pvalue <- paste0('(P=', adonis(as.dist(dist.temp) ~ grp)$aov.tab[1, 6], ')')
    } else {
      pvalue <- ''
    }

    colnames(y) <- c("PC1", "PC2")
    xlim <- c(-max(range(abs(y[, 1]))) * 1.25, max(range(abs(y[, 1]))) * 1.25)
    ylim <- c(-max(range(abs(y[, 2]))) * 1.25, max(range(abs(y[, 2]))) * 1.25)
    plot(y[, 1], y[, 2], type='n', xlim=xlim, ylim=ylim,
         xlab=xlab, ylab=ylab)
    #		points(y[, 1], y[, 2], bg='yellow', col=darkcols[grp],
    #				type='p', pch = pchs[grp], cex=cex.pt)
    col1 <- lightcols[1:nlevels(grp)]
    col2 <- darkcols[1:nlevels(grp)]
    col3 <- darkcols[grp]
    if (!is.null(emp.lev)) {

      col1 <- col2
      col1[which(levels(grp) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
      col1[which(levels(grp) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
      col2[which(levels(grp) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
      col2[which(levels(grp) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
      col3[grp !=  emp.lev] <- rgb(0.5, 0.5, 0.5, 0.5)
      col3[grp ==  emp.lev] <- rgb(1.0, 0, 0, 1.0)

    }
    if (ellipse == T) {
      s.class(y,
              fac = grp,
              cstar = cstar,
              clab = 0,
              cpoint = 0,
              axesell = F,
              col = col1,
              grid = TRUE,
              add.plot=T,
              ...
      )
    }

    points(y[, 1], y[, 2], bg=col3, col='black',
           type='p', pch = pchs[strata], cex=cex.pt)
    if (indiv.lab) {
      if (is.null(emp.lev)) {
        lab.temp <- rownames(dist.temp)
        text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
      } else {
        lab.temp <- rownames(dist.temp)
        lab.temp[grp !=  emp.lev] <- ''
        text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
      }

    }
    s.class(y,
            fac = grp,
            cstar =0,
            cellipse = 0,
            clab = clab,
            cpoint = 0,
            axesell = F,
            col = col2,
            grid = TRUE,
            add.plot=T
    )
    if (!is.null(strata0)) {
      legend('topright', legend=(levels(strata)), pch=pchs[1:nlevels(strata)])
    }

    title(main=paste(dist.name, "distance"), sub=paste0(sub, pvalue))
    #		text(-0.25, 0.3, "PERMANOVA p=0.016")
  }

  if (pdf) {
    dev.off()
  }
}


# New: 2018_07_16  PCs trend
generate_PC_scatterplot <- function (data.obj, dist.obj=NULL, grp.name, adj.name=NULL, subject=NULL,
                                     dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), nPC=2, pca.method='cmd',
                                     strata=NULL, strata2=NULL, combine.strata=FALSE, combine.strata2=FALSE, smooth.method='loess',
                                     pt.shape=16, pt.alpha=0.5, pt.size=2, subject.pt = FALSE,
                                     ann='', hei=NULL, wid=NULL, pdf=TRUE, gg.cmd=NULL) {

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

    if (!is.null(strata) & combine.strata == FALSE) {
      wid <- 4 * nlevels(factor(df[, strata]))
    }
    if (!is.null(strata2) & combine.strata2 == FALSE) {
      hei <- 2.5 * nlevels(factor(df[, strata2]))
    }
  }

  if (pdf == TRUE) {
    pdf(paste0('Beta_diversity_PC_scatterplot_', ann, '.pdf'), height=hei, width=wid)
  }

  measures <- NULL
  x <- NULL
  for (dist.name in dist.names) {
    dist.temp <- dist.obj[[dist.name]]
    if (!is.null(adj.name)) {
      adj <- as.data.frame(df[, adj.name])
      obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
      dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
      dist.temp <- dist(dat2)
    }
    if (pca.method == 'cmd') {
      obj <- cmdscale(as.dist(dist.temp), k=nPC, eig=T)
      pve <- round(obj$eig[1:nPC]/sum(abs(obj$eig))*100, 1)
      x <- cbind(x, obj$points)
      measures <- c(measures, paste0(dist.name, ' - MDS PC', 1:nPC, '(', pve, '%)'))
    }

    if (pca.method == 'nmds') {
      obj <- metaMDS(as.dist(dist.temp), k=nPC)
      x <- cbind(x, obj$points)
      measures <- c(measures, paste0(dist.name, ' - NMDS PC', 1:nPC))
    }
  }
  colnames(x) <- measures


  if (is.null(strata)) {
    for (measure in measures) {
      cat(measure, '\n')
      xx <- x[, measure]
      df2 <- data.frame(Value=xx, Group=grp, ID=ID)
      obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
        geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
        labs(y='Value', x=grp.name, title=measure) +
        theme(legend.position="none")
      if (!is.null(gg.cmd)) {
        obj <- obj + eval(parse(text=gg.cmd))
      }

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
    }
  } else {
    for (measure in measures) {
      cat(measure, '\n')

      xx <- x[, measure]
      grp2 <- factor(df[, strata])
      if (!is.null(strata2)) {
        grp3 <- factor(df[, strata2])
      } else {
        grp3 <- grp2
      }
      grp4 <- factor(paste(grp2, grp3))

      # VERY bad names with only difference in capital letters, variable names are also confusing
      # Rev: 2018_07_16

      df2 <- data.frame(Value=xx, Group=grp, Strata=grp2, Strata2=grp3, Strata3=grp4, ID=ID)

      if (is.null(strata2)) {
        if (combine.strata == TRUE) {
          obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata)) +
            #		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
            geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
            labs(y='Value', x=grp.name, title=measure)
        } else {
          obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
            #		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
            geom_smooth(method=smooth.method,  aes(group=NA),  size=0.75) +
            labs(y='Value', x=grp.name, title=measure) +
            facet_grid(. ~ Strata)
        }

      } else {
        if (combine.strata == TRUE) {
          if (combine.strata2 == TRUE) {
            obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata3)) +
              #			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
              geom_smooth(method=smooth.method, aes(group=Strata3), size=0.75) +
              labs(y='Value', x=grp.name, title=measure)
          } else {
            obj <- ggplot(df2, aes(x=Group, y=Value, group=ID, col=Strata)) +
              #			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
              geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
              labs(y='Value', x=grp.name, title=measure) +
              facet_grid(Strata2 ~ .)
          }
        } else {
          if (combine.strata2 == TRUE) {
            warning('Currently, combine.strata2 will be forced to be FALSE if combine.strata=FALSE!\n')
          }
          obj <- ggplot(df2, aes(x=Group, y=Value, group=ID)) +
            #	geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
            geom_smooth(method=smooth.method,  aes(group=NA), size=0.75) +
            labs(y='Value', x=grp.name, title=measure) +
            facet_grid(Strata2 ~ Strata)
        }

      }

      obj <- obj + theme(legend.title=element_blank())

      if (!is.null(gg.cmd)) {
        obj <- obj + eval(parse(text=gg.cmd))
      }

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
    }

  }
  if (pdf == TRUE) {
    dev.off()
  }


}


# New: 2018_07_16  grp.name codes the day
generate_ordination_trajectory <- function (data.obj, dist.obj=NULL, day.name=NULL, adj.name=NULL, subject=NULL,
                                            dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),  pca.method='cmd',
                                            strata=NULL, strata2=NULL, combine.strata=FALSE, combine.strata2=FALSE, smooth.method='auto',
                                            pt.shape=16, pt.alpha=0.75, pt.size=2, subject.pt = TRUE,
                                            ann='', hei=NULL, wid=NULL, pdf=TRUE, gg.cmd=NULL) {

  nPC <- 2
  df <- data.obj$meta.dat
  ind <- order(df[, day.name])

  df <- df[ind, ]
  dist.obj <- lapply(dist.obj, function (x) x[ind, ind])

  day <- df[, day.name]


  if(is.null(subject)) {
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
    pdf(paste0('Beta_diversity_Subject_Trajectory_', ann, '.pdf'), height=hei, width=wid)
  }

  x1.measures <- NULL
  x2.measures <- NULL
  x1 <- NULL
  x2 <- NULL
  for (dist.name in dist.names) {
    dist.temp <- dist.obj[[dist.name]]
    if (!is.null(adj.name)) {
      adj <- as.data.frame(df[, adj.name])
      obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
      dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
      dist.temp <- dist(dat2)
    }
    if (pca.method == 'cmd') {
      obj <- cmdscale(as.dist(dist.temp), k=nPC, eig=T)
      pve <- round(obj$eig[1:nPC]/sum(abs(obj$eig))*100, 1)
      x1 <- cbind(x1, obj$points[, 1])
      x2 <- cbind(x2, obj$points[, 2])
      x1.measures <- c(x1.measures, paste0('MDS PC', 1, '(', pve[1], '%)'))
      x2.measures <- c(x2.measures, paste0('MDS PC', 2, '(', pve[2], '%)'))
    }

    if (pca.method == 'nmds') {
      obj <- metaMDS(as.dist(dist.temp), k=nPC)
      x1 <- cbind(x1, obj$points[, 1])
      x2 <- cbind(x2, obj$points[, 2])
      x1.measures <- c(x1.measures, paste0('MDS PC', 1))
      x2.measures <- c(x2.measures, paste0('MDS PC', 2))
    }
  }
  colnames(x1) <- dist.names
  colnames(x2) <- dist.names
  names(x1.measures) <- dist.names
  names(x2.measures) <- dist.names
  measures <- dist.names


  if (is.null(strata)) {
    for (measure in measures) {
      cat(measure, '\n')
      xx1 <- x1[, measure]
      xx2 <- x2[, measure]
      df2 <- data.frame(Day=factor(day), ID=ID, xx1=xx1, xx2=xx2)
      obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day)) +
        #		geom_smooth(method=smooth.method, aes(group=NA), size=0.75) +
        labs(y=x2.measures[measure], x=x1.measures[measure], title=measure)


      if (subject.pt == FALSE) {
        obj <- obj + geom_path(size=0.25, alpha=0.5, colour='black')
      }  else {
        obj <- obj + geom_path(size=0.25, alpha=0.5, colour='black') + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha)
      }

      if (!is.null(gg.cmd)) {
        obj <- obj + eval(parse(text=gg.cmd))
      }

      obj <- obj + theme(legend.title=element_blank())

      print(obj)
    }
  } else {
    for (measure in measures) {
      cat(measure, '\n')

      xx1 <- x1[, measure]
      xx2 <- x2[, measure]
      grp2 <- factor(df[, strata])
      if (!is.null(strata2)) {
        grp3 <- factor(df[, strata2])
      } else {
        grp3 <- grp2
      }
      grp4 <- factor(paste(grp2, grp3))

      # VERY bad names with only difference in capital letters, variable names are also confusing
      # Rev: 2018_07_16

      df2 <- data.frame(xx1=xx1, xx2=xx2, Strata=grp2, Strata2=grp3, Strata3=grp4, ID=ID, Day=factor(day))

      if (is.null(strata2)) {
        if (combine.strata == TRUE) {
          obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day, linetype=Strata)) +
            #		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
            #	geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
            labs(y=x2.measures[measure], x=x1.measures[measure], title=measure)

          if (subject.pt == FALSE) {
            obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata))
          }  else {
            obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata)) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha)
          }

        } else {
          obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day)) +
            #		geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
            # geom_smooth(method=smooth.method,  aes(group=NA),  size=0.75) +
            labs(y=x2.measures[measure], x=x1.measures[measure], title=measure) +
            facet_grid(. ~ Strata)
          if (subject.pt == FALSE) {
            obj <- obj + geom_path(size=0.25, alpha=0.5, col='black')
          }  else {
            obj <- obj + geom_path(size=0.25, alpha=0.5, col='black') + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha)
          }

        }

      } else {
        if (combine.strata == TRUE) {
          if (combine.strata2 == TRUE) {
            obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day, linetype=Strata3)) +
              #			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
              # geom_smooth(method=smooth.method, aes(group=Strata3), size=0.75) +
              labs(y=x2.measures[measure], x=x1.measures[measure], title=measure)
            if (subject.pt == FALSE) {
              obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata3))
            }  else {
              obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata3)) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha)
            }
          } else {
            obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day, linetype=Strata)) +
              #			geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
              # geom_smooth(method=smooth.method, aes(group=Strata), size=0.75) +
              labs(y=x2.measures[measure], x=x1.measures[measure], title=measure) +
              facet_grid(Strata2 ~ .)
            if (subject.pt == FALSE) {
              obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata))
            }  else {
              obj <- obj + geom_path(size=0.25, alpha=0.75, aes(col=Strata)) + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha)
            }

          }
        } else {
          if (combine.strata2 == TRUE) {
            warning('Currently, combine.strata2 will be forced to be FALSE if combine.strata=FALSE!\n')
          }
          obj <- ggplot(df2, aes(x=xx1, y=xx2, group=ID, col=Day)) +
            #	geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha) +
            #	geom_smooth(method=smooth.method,  aes(group=NA), size=0.75) +
            labs(y=x2.measures[measure], x=x1.measures[measure], title=measure) +
            facet_grid(Strata2 ~ Strata)
          if (subject.pt == FALSE) {
            obj <- obj + geom_path(size=0.25, alpha=0.5, colour='black')
          }  else {
            obj <- obj + geom_path(size=0.25, alpha=0.5, colour='black') + geom_point(size=pt.size, shape=pt.shape, alpha=pt.alpha)
          }
        }

      }

      obj <- obj + theme(legend.title=element_blank())

      if (!is.null(gg.cmd)) {
        obj <- obj + eval(parse(text=gg.cmd))
      }

      print(obj)
    }

  }
  if (pdf == TRUE) {
    dev.off()
  }


}

# New: 2017_02_21  Separate PC plots for strata2
# Rev: 2017_05_04 Add p value
# Rev: 2017_06_01 Different pages for different distances
generate_ordination_separate <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                          grp.name, adj.name=NULL, emp.lev=NULL, indiv.lab=FALSE, indiv.lab.cex=0.5, is.lightcolor=TRUE,
                                          strata=NULL,  separate=NULL,  layout=NULL, pc.separate, pca.method='cmd', ann=NULL, sub=NULL,
                                          clab=1.0, cex.pt=1.25, ellipse=T, cstar= 1, wid=NULL, hei=NULL, pdf=TRUE, is.pvalue = FALSE, ...) {
  # Implment strata
  # To be completed, add continuous case
  strata0 <- strata
  if (is.null(separate)){
    stop('The augment separate is required!\n')
  } else {
    data.obj$meta.dat[, separate] <- factor(data.obj$meta.dat[, separate])
    sep.levels <- levels(data.obj$meta.dat[, separate])
  }

  if (!is.null(layout)) {
    if (is.null(hei) | is.null(wid)) {
      hei <- 3 * layout[1]
      wid <- 3 * layout[2]
    }

  } else {
    layout <- c(1, length(sep.levels))
    if (is.null(hei) | is.null(wid)) {
      hei <- 3
      wid <- 3 * length(sep.levels)
    }
  }

  if (pdf) {

    if (is.null(ann)) {
      pdf(paste0('Beta_diversity_ordination_', pca.method, '_', grp.name, '_separate.pdf'), width=wid, height=hei)
    } else {
      pdf(paste0('Beta_diversity_ordination_', pca.method, '_', ann, '_separate.pdf'), width=wid, height=hei)
    }
  }

  par(mfrow=layout)

  df <- data.obj$meta.dat
  grp <- factor(df[, grp.name])

  if (!is.null(emp.lev)) {
    grp <- factor(grp, levels=c(setdiff(levels(grp), emp.lev), emp.lev))
  }

  if (!is.null(strata)) {
    strata <- factor(df[, strata])
  } else {
    strata <- factor(grp)
  }

  # To be revised
  darkcols <- hue_pal(l=40)(nlevels(grp))
  if (is.lightcolor) {
    lightcols <- hue_pal(c=45, l=80)(nlevels(grp))
  } else {
    lightcols <- darkcols
  }

  pchs <- rep(c(21, 22, 23, 24, 25), ceiling(nlevels(strata) / 5))[1:nlevels(strata)]

  for (dist.name in dist.names) {
    dist.temp <- dist.obj[[dist.name]]


    for (sep.level in sep.levels) {
      ind <- df[, separate] == sep.level
      df2 <- df[ind, ]
      grp2 <- grp[ind]
      strata2 <- strata[ind]
      dist.temp2 <- dist.temp[ind, ind]

      if (!is.null(adj.name)) {
        adj <- as.data.frame(df2[, adj.name])
        obj <- cmdscale(as.dist(dist.temp2), k=ncol(dist.temp2)-1)
        dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
        dist.temp2 <- dist(dat2)
      }
      if (pca.method == 'cmd') {
        obj <- cmdscale(as.dist(dist.temp2), k=2, eig=T)
        pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
        y <- cbind(obj$points[, 1], obj$points[, 2])

        xlab <- paste0('PC1(', pve[1], '%)')
        ylab <- paste0('PC2(', pve[2], '%)')
      }

      if (pca.method == 'nmds') {
        obj <- metaMDS(as.dist(dist.temp2), k=2)
        y <- cbind(obj$points[, 1], obj$points[, 2])
        xlab <- 'NMDS1'
        ylab <- 'NMDS2'
      }

      if (pca.method == 'pls') {
        require(mixOmics)
        # Test
        obj <- cmdscale(as.dist(dist.temp2), k=ncol(dist.temp2)-1)
        obj <- plsda(obj, grp2, ncomp=2)
        y <- cbind(obj$variates$X[, 1], obj$variates$X[, 2])
        xlab <- 'PLS1'
        ylab <- 'PLS2'
      }

      if (is.pvalue == TRUE & nlevels(factor(grp2)) > 1) {
        pvalue <- paste0('(P=', adonis(as.dist(dist.temp2) ~ grp2)$aov.tab[1, 6], ')')
      } else {
        pvalue <- ''
      }
      colnames(y) <- c("PC1", "PC2")
      xlim <- c(-max(range(abs(y[, 1]))) * 1.25, max(range(abs(y[, 1]))) * 1.25)
      ylim <- c(-max(range(abs(y[, 2]))) * 1.25, max(range(abs(y[, 2]))) * 1.25)
      plot(y[, 1], y[, 2], type='n', xlim=xlim, ylim=ylim,
           xlab=xlab, ylab=ylab)
      #		points(y[, 1], y[, 2], bg='yellow', col=darkcols[grp2],
      #				type='p', pch = pchs[grp2], cex=cex.pt)
      col1 <- lightcols[1:nlevels(grp2)]
      col2 <- darkcols[1:nlevels(grp2)]
      col3 <- darkcols[grp2]
      if (!is.null(emp.lev)) {

        col1 <- col2
        col1[which(levels(grp2) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
        col1[which(levels(grp2) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
        col2[which(levels(grp2) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
        col2[which(levels(grp2) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
        col3[grp2 !=  emp.lev] <- rgb(0.5, 0.5, 0.5, 0.5)
        col3[grp2 ==  emp.lev] <- rgb(1.0, 0, 0, 1.0)

      }
      if (ellipse == T) {
        s.class(y,
                fac = grp2,
                cstar = cstar,
                clab = 0,
                cpoint = 0,
                axesell = F,
                col = col1,
                grid = TRUE,
                add.plot=T,
                ...
        )
      }

      points(y[, 1], y[, 2], bg=col3, col='black',
             type='p', pch = pchs[strata2], cex=cex.pt)
      if (indiv.lab) {
        if (is.null(emp.lev)) {
          lab.temp <- rownames(dist.temp2)
          text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
        } else {
          lab.temp <- rownames(dist.temp2)
          lab.temp[grp2 !=  emp.lev] <- ''
          text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
        }

      }
      s.class(y,
              fac = grp2,
              cstar =0,
              cellipse = 0,
              clab = clab,
              cpoint = 0,
              axesell = F,
              col = col2,
              grid = TRUE,
              add.plot=T
      )
      if (!is.null(strata0)) {
        legend('topright', legend=(levels(strata2)), pch=pchs[1:nlevels(strata2)])
      }

      title(main=sep.level, sub=paste(dist.name, "distance", pvalue))

    }


    #		text(-0.25, 0.3, "PERMANOVA p=0.016")
  }

  if (pdf) {
    dev.off()
  }
}


# Rev:2017_04_19 Add bootstrap standard error

generate_distance_barplot <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                       grp.name, strata=NULL, within=T, between=T, bt.no = 100, ann='') {
  strata.name <- strata
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  grp.levels <- levels(grp)
  grp.nlevels <- nlevels(grp)
  grp.btws <- outer(grp.levels, grp.levels, paste, sep="#")
  grp.btws <- grp.btws[lower.tri(grp.btws)]

  if (!is.null(strata.name)) {
    strata <- df[, strata.name]
  } else {
    strata <- factor(rep(1, nrow(df))) # pseudo strata
  }
  res.df <- NULL
  for (dist.name in dist.names) {
    for (stratum in levels(strata)) {
      ind <- strata %in% stratum
      dist.sub <- dist.obj[[dist.name]][ind, ind]
      df2 <- df[ind, , drop=FALSE]
      if (between) {
        for (grp.btw in grp.btws) {
          ind1 <- which(df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[1])
          ind2 <- which(df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[2])
          temp <- as.vector(dist.sub[ind1, ind2])
          sem <- sd(sapply(1:bt.no, function (i) {
            mean(dist.sub[as.numeric(sample(paste(ind1), repl=TRUE)),
                          as.numeric(sample(paste(ind2), repl=TRUE))])
          }))
          res.df <- rbind(res.df, data.frame(DistanceMetric=dist.name, Strata=stratum, Compare='Between',
                                             DistanceType=grp.btw, Distance=mean(temp), sd=sem))
        }
      }

      if (within) {
        for (grp.wth in grp.levels) {
          ind1 <- which(df2[, grp.name] == grp.wth)
          temp <- dist.sub[ind1, ind1]
          temp <- temp[lower.tri(temp)]
          sem <- sd(sapply(1:bt.no, function (i) {
            ind2 <- as.numeric(sample(paste(ind1), repl = TRUE))
            temp <- dist.sub[ind2, ind2]
            temp <- temp[lower.tri(temp)]
            mean(temp)
          }))
          res.df <- rbind(res.df, data.frame(DistanceMetric=dist.name, Strata=stratum, Compare='Within',
                                             DistanceType=grp.wth, Distance=mean(temp), sd=sem))
        }
      }
    }

  }
  if (between & within) {
    res.df$DistanceType <- factor(res.df$DistanceType, levels=c(grp.btws, grp.levels))
    levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'), paste('Within\n', grp.levels))
  } else {
    if (between) {
      res.df$DistanceType <- factor(res.df$DistanceType)
      levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'))
    } else {
      res.df$DistanceType <- factor(res.df$DistanceType)
      levels(res.df$DistanceType) <- c(paste('Within\n', grp.levels))
    }
  }

  if (is.null(strata.name)) {
    pdf(paste0("Beta_diversity_btw", between, "_wth", within, "no_strata_barplot_", ann, ".pdf"),
        width=5, height=5)

    limits <- aes(ymax = Distance + sd, ymin=Distance - sd)
    dodge <- position_dodge(width=0.9)
    for (dist.name in dist.names) {
      cat(dist.name, "...\n")
      temp <- res.df[res.df$DistanceMetric == dist.name, ]
      obj1 <- ggplot(temp, aes(x=DistanceType, y=Distance, fill=DistanceType)) +
        geom_bar(position=dodge, stat="identity", width=0.75) +
        geom_bar(position=dodge, stat="identity", width=0.75, colour="black", show_guide=FALSE, size=0.25) +
        geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
        labs(y=paste(dist.name, "Distance"), x='') +
        theme(legend.position="none") +
        theme(axis.text.x=element_text(angle=90, hjust=1))

      print(obj1)

    }
    dev.off()
  } else {
    pdf(paste0("Beta_diversity_btw", between, "_wth", within, "_strata", strata.name, "barplot_", ann, ".pdf"),
        width=2.5*(nlevels(strata)-1) + 5, height=5)

    limits <- aes(ymax = Distance + sd, ymin=Distance - sd)
    dodge <- position_dodge(width=0.9)
    for (dist.name in dist.names) {
      cat(dist.name, "...\n")
      temp <- res.df[res.df$DistanceMetric == dist.name, ]
      obj1 <- ggplot(temp, aes(x=Strata, y=Distance, fill=DistanceType)) +
        geom_bar(position=dodge, stat="identity") +
        geom_bar(position=dodge, stat="identity", colour="black", show_guide=FALSE, size=0.25) +
        geom_errorbar(limits, position=dodge, size=0.25, width=0.5) +
        labs(y=paste(dist.name, "Distance"), x=strata.name) +
        theme(axis.text.x=element_text(angle=90, hjust=1))
      print(obj1)

    }
    dev.off()

  }

}

generate_distance_boxplot <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                       grp.name, strata=NULL, within=F, between=T, ann='') {
  strata.name <- strata
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  grp.levels <- levels(grp)
  grp.nlevels <- nlevels(grp)
  grp.btws <- outer(grp.levels, grp.levels, paste, sep="#")
  grp.btws <- grp.btws[lower.tri(grp.btws)]

  if (!is.null(strata.name)) {
    strata <- df[, strata.name]
  } else {
    strata <- factor(rep(1, nrow(df))) # pseudo strata
  }
  res.df <- NULL
  for (dist.name in dist.names) {
    for (stratum in levels(strata)) {
      ind <- strata %in% stratum
      dist.sub <- dist.obj[[dist.name]][ind, ind]
      df2 <- df[ind, , drop=FALSE]
      if (between) {
        for (grp.btw in grp.btws) {
          ind1 <- df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[1]
          ind2 <- df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[2]
          temp <- as.vector(dist.sub[ind1, ind2])
          n <- length(temp)
          res.df <- rbind(res.df, data.frame(DistanceMetric=rep(dist.name, n), Strata=rep(stratum, n), Compare=rep('Between', n),
                                             DistanceType=rep(grp.btw, n), Distance=temp))
        }
      }

      if (within) {
        for (grp.wth in grp.levels) {
          ind1 <- df2[, grp.name] == grp.wth
          temp <- dist.sub[ind1, ind1]
          temp <- temp[lower.tri(temp)]
          n <- length(temp)
          res.df <- rbind(res.df, data.frame(DistanceMetric=rep(dist.name, n), Strata=rep(stratum, n), Compare=rep('Within', n),
                                             DistanceType=rep(grp.wth, n), Distance=temp))
        }
      }
    }

  }
  if (between & within) {
    res.df$DistanceType <- factor(res.df$DistanceType, levels=c(grp.btws, grp.levels))
    levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'), paste('Within\n', grp.levels))
  } else {
    if (between) {
      res.df$DistanceType <- factor(res.df$DistanceType)
      levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'))
    } else {
      res.df$DistanceType <- factor(res.df$DistanceType)
      levels(res.df$DistanceType) <- c(paste('Within\n', grp.levels))
    }
  }

  if (is.null(strata.name)) {
    pdf(paste0("Beta_diversity_btw", between, "_wth", within, "no_strata_boxplot_", ann, ".pdf"),
        width=5, height=5)

    for (dist.name in dist.names) {
      cat(dist.name, "...\n")
      temp <- res.df[res.df$DistanceMetric == dist.name, ]

      dodge <- position_dodge(width=0.95)
      obj1 <- ggplot(temp, aes(x=DistanceType, y=Distance, col=DistanceType)) +
        geom_boxplot(position=dodge, outlier.colour = NA) +
        #					geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
        labs(y=paste(dist.name, "Distance"), x='') +
        theme(legend.position="none")

      print(obj1)

    }
    dev.off()
  } else {
    pdf(paste0("Beta_diversity_btw", between, "_wth", within, "_strata", strata.name, "boxplot_", ann, ".pdf"),
        width=2.5*(nlevels(strata)-1) + 5, height=5)

    for (dist.name in dist.names) {
      cat(dist.name, "...\n")
      temp <- res.df[res.df$DistanceMetric == dist.name, ]
      dodge <- position_dodge(width=0.95)
      obj1 <- ggplot(temp, aes(x=Strata, y=Distance, col=DistanceType)) +
        geom_boxplot(position=dodge, outlier.colour = NA) +
        #			geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
        labs(y=paste(dist.name, "Distance"), x=strata.name)
      print(obj1)

    }
    dev.off()

  }

}

# Rev:2016_11_25
generate_clustering <- function(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), meta.info, cluster.method='average',
                                is.labRow=F, cex.lab=NULL, ann="", wid=10, hei=6, colFnsF=NULL, colFnsC=NULL) {
  if (is.null(colFnsC)) {
    colFnsC <- c(colorRampPalette(c('black', 'green')), colorRampPalette(c('black', 'blue')), colorRampPalette(c('black', 'red')))
  }
  if (is.null(colFnsF)) {
    rainbow3 <- function(x) {
      rainbow(x + 2)[1:x]
    }
    jet3 <- function(x) {
      jet(x + 2)[1:x]
    }
    colFnsF <- c(rainbow3, jet3)
  }

  pdf(paste0("Beta_diversity_Hierachical_clustering_", ann, ".pdf"), width=wid, height=hei)
  for (dist.name in dist.names) {

    df <- data.obj$meta.dat
    dist.sub <- dist.obj[[dist.name]]
    dend <- hclust(as.dist(dist.sub), cluster.method)


    key.list <- list()
    mat <- NULL
    for (keyID in meta.info) {
      x <- df[, keyID]
      i <- 0
      j <- 0
      if (is.factor(x)) {
        key.list[[keyID]] <- list(breaks=levels(x), colors=(colFnsF[i+1][[1]])(nlevels(x)), base=NA, col.na=NA, right=F, include.lowest=F)
        i <- (i + 1) %% length(colFnsF)
        mat <- cbind(mat, key.list[[keyID]]$colors[x])
      } else {
        key.list[[keyID]] <- makecmap(x, n=5, colFn=colFnsC[j+1][[1]])
        j <- (j + 1) %% length(colFnsC)
        mat <- cbind(mat, cmap(x, key.list[[keyID]]))
      }
    }
    colnames(mat) <- meta.info

    par(oma = c(1, 2, 1, 2))
    par(mar = c(5,4,4,3)+0.1)  # make space for color keys
    if (is.labRow == F) {
      labRow <- ""
    } else {
      labRow <- rownames(df)
    }
    if (is.null(cex.lab)) {
      cex.lab <- 25 / nrow(df)
    }

    dendromat(dend, mat, labRow=labRow,
              ylab = 'Distance', main = paste(dist.name, "distance"), cex.lab=cex.lab)

    par(oma=c(0, 0, 0, 0))
    y.cord <- (1/length(meta.info)) * (0:(length(meta.info) - 1))
    k <- 1
    for (keyID in meta.info) {
      x <- df[, keyID]
      if (is.factor(x)) {
        vkey2(key.list[[keyID]], keyID, y=y.cord[k], stretch=1.2)
      } else {
        vkey(key.list[[keyID]], keyID, y=y.cord[k], stretch=1.2)
      }
      k <- k + 1
    }
  }
  dev.off()
}

perform_permanova_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                    PermanovaG.dist=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                    formula=NULL,  grp.name=NULL, adj.name=NULL, pairwise=F, block.perm=F, strata=NULL, ann='', ...) {
  # PermanovaG not implemented for block permutation
  result <- list()

  df <- data.obj$meta.dat
  if (!is.null(strata)) {
    if (is.character(strata)) {
      strata <- factor(df[, strata])
    }
  }
  if (!is.null(formula)) {
    ind <- apply(df[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
    df <- df[ind, ]
    if (!is.null(strata)) {
      strata <- factor(strata[ind])
    }
    sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '.txt'))
    date()
    cat('\nPERMANOVA test: \n')
    permanova.obj <- list()
    for (dist.name in dist.names) {
      cat(dist.name, " distance: \n")
      dist.mat <- as.dist(dist.obj[[dist.name]][ind, ind])
      if (block.perm == F) {
        obj <- adonis(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
      } else {
        obj <- adonis2(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
      }

      prmatrix(obj$aov.tab)
      permanova.obj[[dist.name]] <- obj$aov.tab
      cat("\n")
    }
    result$permanova.obj <- permanova.obj
    permanovaG.obj <- NULL
    if (block.perm == F & !is.null(PermanovaG.dist)) {
      cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
      response <- array(NA, c(sum(ind), sum(ind), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
      for (dist.name in PermanovaG.dist) {
        response[, , dist.name] <- dist.obj[[dist.name]][ind, ind]
      }
      obj <- PermanovaG2(as.formula(paste("response", formula)), df,  strata=strata, ...)
      prmatrix(obj$aov.tab)
      permanovaG.obj <- obj$aov.tab
      cat("\n")
      result$permanovaG.obj <- permanovaG.obj
    }
    cat("\n")
    sink()
  } else {
    if (pairwise == F) {
      if (is.null(adj.name)) {
        formula <- paste('~', grp.name)
      } else {
        formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)
      }

      ind <- apply(df[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
      df <- df[ind, ]
      if (!is.null(strata)) {
        strata <- factor(strata[ind])
      }
      sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '.txt'))
      date()
      cat('\nPERMANOVA test: \n')
      permanova.obj <- list()
      for (dist.name in dist.names) {
        cat(dist.name, " distance: \n")
        dist.mat <- as.dist(dist.obj[[dist.name]][ind, ind])
        if (block.perm == F) {
          obj <- adonis(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
        } else {
          obj <- adonis2(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
        }
        prmatrix(obj$aov.tab)
        permanova.obj[[dist.name]] <- obj$aov.tab
        cat("\n")
      }
      result$permanova.obj <- permanova.obj
      permanovaG.obj <- NULL
      if (block.perm == F & !is.null(PermanovaG.dist)) {
        cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
        response <- array(NA, c(sum(ind), sum(ind), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
        for (dist.name in PermanovaG.dist) {
          response[, , dist.name] <- dist.obj[[dist.name]][ind, ind]
        }
        obj <- PermanovaG2(as.formula(paste("response", formula)), df,  strata=strata, ...)
        prmatrix(obj$aov.tab)
        permanovaG.obj <- obj$aov.tab
        cat("\n")
        result$permanovaG.obj <- permanovaG.obj
      }

      cat("\n")
      sink()

    } else {
      sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwise.txt'))
      date()
      cat('\nPairwise PERMANOVA test: \n')
      grp <- factor(df[, grp.name])
      grp.levels <- levels(grp)
      grp.nlevels <- nlevels(grp)
      pmat.all <- NULL
      rmat.all <- NULL
      pmat.G <- matrix(NA, grp.nlevels, grp.nlevels)
      colnames(pmat.G) <- rownames(pmat.G) <- grp.levels
      for (dist.name in dist.names) {
        cat(dist.name, " distance: \n")
        pmat <- matrix(NA, grp.nlevels, grp.nlevels)
        colnames(pmat) <- rownames(pmat) <- grp.levels
        rmat <- matrix(NA, grp.nlevels, grp.nlevels)
        colnames(rmat) <- rownames(rmat) <- grp.levels
        for (i in 1:(grp.nlevels-1)) {
          grp.level1 <- grp.levels[i]
          for (j in (i+1):grp.nlevels) {

            grp.level2 <- grp.levels[j]
            cat(grp.level1, ' vs ', grp.level2, '\n')
            ind <- grp %in% c(grp.level1, grp.level2)
            df2 <- subset(df, ind)
            df2[, grp.name] <- factor(df2[, grp.name])
            dist.mat <- dist.obj[[dist.name]][ind, ind]

            if (!is.null(strata)) {
              strata2 <- factor(strata[ind])
            } else {
              strata2 <- NULL
            }

            if (is.null(adj.name)) {
              formula <- paste('~', grp.name)
            } else {
              formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)
            }

            ind2 <- apply(df2[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
            df2 <- df2[ind2, ]
            dist.mat2 <- as.dist(dist.mat[ind2, ind2])

            if (!is.null(strata2)) {
              strata2 <- factor(strata2[ind2])
            }

            if (block.perm == F) {
              obj <- adonis(as.formula(paste("dist.mat2", formula)), df2,  strata=strata2, ...)
            } else {
              obj <- adonis2(as.formula(paste("dist.mat2", formula)), df2,  strata=strata2,  ...)
            }
            prmatrix(obj$aov.tab)
            cat("\n")

            if (block.perm == F) {
              pmat[i, j] <- pmat[j, i] <- obj$aov.tab[length(adj.name)+1, 6]
              rmat[i, j] <- rmat[j, i] <- obj$aov.tab[length(adj.name)+1, 5]
            } else {
              pmat[i, j] <- pmat[j, i] <- obj$aov.tab[1, 6]
              rmat[i, j] <- rmat[j, i] <- obj$aov.tab[1, 5]
            }

            # PERMANOVA G after last distance
            if (block.perm == F & !is.null(PermanovaG.dist)) {
              if (dist.name == dist.names[length(dist.names)]) {
                cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
                response <- array(NA, c(sum(ind2), sum(ind2), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
                for (dist.name in PermanovaG.dist) {
                  # Rev: 2019_05_20 Fixed an error
                  response[, , dist.name] <- dist.obj[[dist.name]][ind, ind][ind2, ind2]
                }
                obj <- PermanovaG2(as.formula(paste("response", formula)), df2,  strata=strata2, ...)
                prmatrix(obj$aov.tab)
                cat("\n")
                pmat.G[i, j] <- pmat.G[j, i] <- obj$aov.tab[length(adj.name)+1, 'omni.p.value']
              }
            }
          }
        }
        cat("\n")
        pmat.all <- rbind(pmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(pmat), rep("", grp.nlevels))
        rmat.all <- rbind(rmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(rmat), rep("", grp.nlevels))
      }
      cat("\n")
      sink()
      write.csv(pmat.all, paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwiseP.csv'))
      write.csv(rmat.all, paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwiseR.csv'))
      write.csv(pmat.G, paste0('Beta_diversity_PERMANOVA_G_test_', ann, '_pairwiseP.csv'))
    }
  }

  return(invisible(result))
}



# Rev: 2016_12_02, MiKRAT for binary result, add out_type='D'
perform_mirkat_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                 grp.name=NULL, adj.name=NULL, pairwise=F,  ann='', ...) {

  # MiRKAT not implemented for correlated data
  df <- data.obj$meta.dat

  if (pairwise == F) {

    ind <- apply(df[, c(grp.name, adj.name), drop=F], 1, function(x) sum(is.na(x))) == 0
    df <- df[ind, ]

    grp <- df[, grp.name]

    if (is.character(grp)) {
      grp <- factor(grp)
    }
    # Rev: 2016_12_02
    if (is.factor(grp)) {
      if (nlevels(grp) > 2) {
        stop('Currently MiRKAT only supports binary outcome!')
      } else {
        grp <- as.numeric(grp) - 1
        out_type <- 'D'
      }
    } else {
      out_type <- 'C'
    }
    if (!is.null(adj.name)) {
      # No intercept
      adj <- model.matrix(~ ., data.frame(df[, adj.name]))
      # Remove collinear terms
      qadj <- qr(adj, tol = 1e-07)
      adj <- adj[, qadj$pivot, drop = FALSE]
      adj <- adj[, 1:qadj$rank, drop = FALSE]
      # Remove intercept
      adj <- adj[, colSums(adj==1) != nrow(adj)]
    } else {
      adj <- NULL
    }

    sink(paste0('Beta_diversity_MiRKAT_test_', ann, '.txt'))
    date()
    cat('\nMiRKAT test combining ', paste(dist.names, collapse=','), '\n')
    Ks <- list()
    for (dist.name in dist.names) {
      Ks[[dist.name]] <- D2K(dist.obj[[dist.name]][ind, ind])
    }
    # Rev: 2016_12_02
    obj <- MiRKAT(grp, X=adj, Ks, out_type=out_type)
    cat('Individual P value: ')
    prmatrix(t(obj$indivP))
    cat('\nOmnibus P value: ')
    cat(obj$omnibus_p)
    cat("\n")
    sink()

  } else {
    sink(paste0('Beta_diversity_MiRKAT_test_', ann, '_pairwise.txt'))
    date()
    cat('\nPairwise MiRKAT test: \n')
    grp <- factor(df[, grp.name])
    grp.levels <- levels(grp)
    grp.nlevels <- nlevels(grp)
    pmat.G <- matrix(NA, grp.nlevels, grp.nlevels)
    colnames(pmat.G) <- rownames(pmat.G) <- grp.levels
    parr <- array(NA, c(grp.nlevels, grp.nlevels, length(dist.names)), dimnames=list(grp.levels, grp.levels, dist.names))

    for (i in 1:(grp.nlevels-1)) {
      grp.level1 <- grp.levels[i]
      for (j in (i+1):grp.nlevels) {

        grp.level2 <- grp.levels[j]
        cat(grp.level1, ' vs ', grp.level2, '\n')
        ind <- grp %in% c(grp.level1, grp.level2)
        df2 <- subset(df, ind)
        df2[, grp.name] <- factor(df2[, grp.name])
        ind2 <- apply(df2[, c(grp.name, adj.name), drop=F], 1, function(x) sum(is.na(x))) == 0
        df2 <- df[ind2, ]

        grp2 <- as.numeric(df2[, grp.name]) - 1
        if (!is.null(adj.name)) {
          # No intercept
          adj <- model.matrix(~ ., data.frame(df2[, adj.name]))
          # Remove collinear terms
          qadj <- qr(adj, tol = 1e-07)
          adj <- adj[, qadj$pivot, drop = FALSE]
          adj <- adj[, 1:qadj$rank, drop = FALSE]
          # Remove intercept
          adj <- adj[, colSums(adj==1) != nrow(adj)]
        } else {
          adj <- NULL
        }

        Ks <- list()
        for (dist.name in dist.names) {
          Ks[[dist.name]] <- D2K(dist.obj[[dist.name]][rownames(df2), rownames(df2)])
        }
        # Rev: 2016_12_02
        obj <- MiRKAT(grp2, X=adj, Ks, out_type='D')
        pmat.G[i, j] <- pmat.G[j, i] <- obj$omnibus_p
        parr[i, j, ] <- parr[j, i, ] <- obj$indivP
        cat('Individual P value: ')
        prmatrix(t(obj$indivP))
        cat('\nOmnibus P value: ')
        cat(obj$omnibus_p)
        cat("\n")

      }
    }
    cat("\n")
    sink()

    pmat.all <- NULL
    for (dist.name in dist.names) {
      pmat <- parr[, , dist.name]
      pmat.all <- rbind(pmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(pmat), rep("", grp.nlevels))
    }

    write.csv(pmat.all, paste0('Beta_diversity_MiRKAT_test_', ann, '_pairwiseP.csv'))
    write.csv(pmat.G, paste0('Beta_diversity_MiRKAT_O_test_', ann, '_pairwiseP.csv'))
  }
}


perform_betadisper_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), grp.name) {
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  sink('Beta_diversity_BETADISPER_test.txt')
  date()
  cat('\nBetadisper test: \n')
  for (dist.name in dist.names) {
    cat(dist.name, " distance: \n")
    dist.mat <- as.dist(dist.obj[[dist.name]])
    obj <- betadisper(dist.mat, grp)
    prmatrix(anova(obj))
    cat("\n")
  }
  cat("\n")
  sink()
}


perform_distance_compare_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                           grp.name, level1, level2, level3, subject=NULL, alternative='greater', nperm=999, seed=123, ann='') {
  sink(paste0('Beta_diversity_test_', level1, '-', level2, '_', alternative, '_', level1, '-', level3, '_', ann, '.txt'))
  set.seed(seed)
  cat('Testing the distance ', level1, '-', level2, ' is ', alternative, ' than ', level1, '-', level3, '\n')
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  IDs <- df[, subject]
  grp.levels <- levels(grp)
  grp.nlevels <- nlevels(grp)
  ind1 <- which(grp == level1)
  ind2 <- which(grp == level2)
  ind3 <- which(grp == level3)
  if (is.null(subject)) {
    ID2 <- NULL
    ID3 <- NULL
  } else {
    ID2 <- as.character(IDs[grp == level2])
    ID3 <- as.character(IDs[grp == level3])
  }
  for (dist.name in dist.names) {
    cat(dist.name, ' Distance:p=')
    pv <- distance_compare_test(dist.obj[[dist.name]], ind1, ind2, ind3, ID2, ID3, alternative, nperm)
    cat(pv, '\n')
  }
  sink()

}
