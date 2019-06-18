perform_sequence_stat_analysis <- function (data.obj, ann='') {
  sink(paste0('Sequence_Analysis_Statistics_', ann, '.txt'))
  otu.tab <- data.obj$otu.tab

  # Sequencing depth
  otu.abund <- rowSums(otu.tab)
  sam.abund <- colSums(otu.tab)
  otu.prev <- rowSums(otu.tab!=0)/ncol(otu.tab)

  otu.abund <- otu.abund[otu.abund >= 1]
  sam.abund <- sam.abund[sam.abund >= 1]
  cat('This data set contains ', length(sam.abund), ' samples after quality controls.\n')
  cat('16S rDNA targeted sequencing yields ', median(sam.abund), 'reads/sample on average (range:', min(sam.abund), '-', max(sam.abund), ').\n')
  cat('Clustering of these 16S sequence tags produces ', sum(otu.abund > 0), ' OTUs at 97% similarity level.\n')

  pdf(paste0('Sequence_Analysis_Statistics_', ann, '.pdf'), height=5, width=5)
  obj <- ggplot2::ggplot(data=data.frame(x=otu.abund), aes(x=x)) + geom_histogram(col='black', fill='gray') + ylab('Frequency') + xlab('Abundance(Total counts)') +
    scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000, 100000))
  print(obj)
  obj <- ggplot2::ggplot(data=data.frame(x=sam.abund), aes(x=x)) + geom_histogram(col='black', fill='gray')  + ylab('Frequency') + xlab('Sequencing depth')
  print(obj)
  obj <- ggplot2::ggplot(data=data.frame(x=otu.prev), aes(x=x))  + ylab('Frequency') + xlab('Prevalence(Occurence frequency)') + geom_histogram(col='black', fill='gray')
  print(obj)
  dev.off()

  phy.abund <- data.obj$abund.list[['Phylum']]
  fam.abund <- data.obj$abund.list[['Family']]
  gen.abund <- data.obj$abund.list[['Genus']]

  phy.prev <- rowSums(phy.abund != 0) / ncol(phy.abund)
  fam.prev <- rowSums(fam.abund != 0) / ncol(phy.abund)
  gen.prev <- rowSums(gen.abund != 0) / ncol(phy.abund)

  phy.abund <- rowMeans(t(t(phy.abund) / sam.abund))
  fam.abund <- rowMeans(t(t(fam.abund) / sam.abund))
  gen.abund <- rowMeans(t(t(gen.abund) / sam.abund))

  cat('These OTUs belong to ', sum(phy.abund > 0), ' phyla,', sum(fam.abund > 0), ' families and ', sum(gen.abund > 0), 'genera.\n\n')

  phy.prev <- sort(phy.prev, decr=T)
  phy.prev <- round(phy.prev[phy.prev >= 0.90] * 100, 2)

  fam.prev <- sort(fam.prev, decr=T)
  fam.prev <- round(fam.prev[fam.prev >= 0.90] * 100, 2)

  gen.prev <- sort(gen.prev, decr=T)
  gen.prev <- round(gen.prev[gen.prev >= 0.90] * 100, 2)

  # Rev: 2017_02_19 ' ' -> '\n'
  cat('\nThe most prevalent phyla are:\n', paste(paste0(names(phy.prev), '(', phy.prev, '%)'), collapse='\n'), '\n')
  cat('\nThe most prevalent families are:\n', paste(paste0(names(fam.prev), '(', fam.prev, '%)'), collapse='\n'), '\n')
  cat('\nand the most prevalent genera are:\n', paste(paste0(names(gen.prev), '(', gen.prev, '%)'), collapse='\n'), '\n\n')

  phy.abund <- sort(phy.abund, decr=T)
  phy.abund <- round(phy.abund[phy.abund >= 0.05] * 100, 2)

  fam.abund <- sort(fam.abund, decr=T)
  fam.abund <- round(fam.abund[fam.abund >= 0.05] * 100, 2)

  gen.abund <- sort(gen.abund, decr=T)
  gen.abund <- round(gen.abund[gen.abund >= 0.05] * 100, 2)

  cat('\nThe most abundant phyla are:\n', paste(paste0(names(phy.abund), '(', phy.abund, '%)'), collapse='\n'), '\n')
  cat('\nThe most abundant families are:\n', paste(paste0(names(fam.abund), '(', fam.abund, '%)'), collapse='\n'), '\n')
  cat('\nand the most abundant genera are:\n', paste(paste0(names(gen.abund), '(', gen.abund, '%)'), collapse='\n'), '\n\n')
  sink()
}
