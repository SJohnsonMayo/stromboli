#' Generate summary string
#'
#' @param data.obj A data object created by load_data
#'
#' @return A character string containing some basic sequencing statistics: the number of samples, median reads per sample, range of reads per sample, number of OTUs/ASVs
#' @export
#'
#' @examples
generate_summary_string <- function(data.obj){
  otu.tab <- data.obj$otu.tab
  # Sequencing depth
  otu.abund <- rowSums(otu.tab)
  sam.abund <- colSums(otu.tab)
  otu.prev <- rowSums(otu.tab != 0)/ncol(otu.tab)
  otu.abund <- otu.abund[otu.abund >= 1]
  sam.abund <- sam.abund[sam.abund >= 1]
  summary_string <- paste0("This data set contains ", length(sam.abund), " samples after quality controls.\n")
  summary_string <- paste0(summary_string, "16S rDNA targeted sequencing yields ", median(sam.abund), "reads/sample on average (range:",
                           min(sam.abund), "-", max(sam.abund), ").\n")
  summary_string <- paste0(summary_string, "Clustering of these 16S sequence tags produces ", sum(otu.abund > 0), " OTUs at 97% similarity level.\n")
  summary_string
}



#' Basic summary plots
#'
#' @param data.obj Data object created by load_data()
#' @param type Type of plot to generate, either "abundance", "depth" or "prevalence"
#'
#' @return A ggplot object based on the type specified by \code{type}
#' @export

#' @examples
#' summary_plot(data.obj, type="abundance")
#'
summary_plot <- function(data.obj, type=""){
  otu.tab <- data.obj$otu.tab
  if(type=="abundance"){
    data <- rowSums(otu.tab)
    obj <- ggplot(data = data.frame(x = data), aes(x = x)) +
      geom_histogram(col = "black", fill = "gray") + ylab("Frequency") +
      xlab("Abundance(Total counts)") +
      scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 1e+05, 1e+05))
  }else if(type=="depth"){
    data <- colSums(otu.tab)
    obj <- ggplot(data = data.frame(x = data), aes(x = x)) +
      geom_histogram(col = "black", fill = "gray") + ylab("Frequency") +
      xlab("Sequencing depth")
  }else if(type=="prevalence"){
    data <- rowSums(otu.tab != 0)/ncol(otu.tab)
    obj <- ggplot(data = data.frame(x = data), aes(x = x)) +
      ylab("Frequency") + xlab("Prevalence(Occurence frequency)") +
      geom_histogram(col = "black", fill = "gray")
  }else{
    return("ERROR: type should be either 'abundance', 'depth', or 'prevalence'")
  }
  obj
}

#' Prevalence and abundance table generation
#'
#' @param data.obj A data object created by load_data()
#'
#' @return A list of tables(phy.prev, fam.prev, gen.prev, phy.abund, fam.abund, and gen.abund) containing the prevalences and abundances for all Phyla, Families, and Genera in \code{data.obj}
#' @export
#'
#' @examples
prev_abund_tables <- function(data.obj){
  sam.abund <- colSums(data.obj$otu.tab)
  phy.abund <- data.obj$abund.list[["Phylum"]]
  fam.abund <- data.obj$abund.list[["Family"]]
  gen.abund <- data.obj$abund.list[["Genus"]]
  phy.prev <- rowSums(phy.abund != 0)/ncol(phy.abund)
  fam.prev <- rowSums(fam.abund != 0)/ncol(phy.abund)
  gen.prev <- rowSums(gen.abund != 0)/ncol(phy.abund)
  phy.abund <- rowMeans(t(t(phy.abund)/sam.abund))
  fam.abund <- rowMeans(t(t(fam.abund)/sam.abund))
  gen.abund <- rowMeans(t(t(gen.abund)/sam.abund))

  list(phy.prev=phy.prev, fam.prev=fam.prev, gen.prev=gen.prev,
       phy.abund=phy.abund, fam.abund=fam.abund, gen.abund=gen.abund)
}

perform_sequence_stat_analysis <- function(data.obj, ann = "") {



  cat("These OTUs belong to ", sum(phy.abund > 0), " phyla,", sum(fam.abund >
                                                                    0), " families and ", sum(gen.abund > 0), "genera.\n\n")
  phy.prev <- sort(phy.prev, decr = T)
  phy.prev <- round(phy.prev[phy.prev >= 0.9] * 100, 2)
  fam.prev <- sort(fam.prev, decr = T)
  fam.prev <- round(fam.prev[fam.prev >= 0.9] * 100, 2)
  gen.prev <- sort(gen.prev, decr = T)
  gen.prev <- round(gen.prev[gen.prev >= 0.9] * 100, 2)
  # Rev: 2017_02_19 ' ' -> '\n'
  cat("\nThe most prevalent phyla are:\n", paste(paste0(names(phy.prev),
                                                        "(", phy.prev, "%)"), collapse = "\n"), "\n")
  cat("\nThe most prevalent families are:\n", paste(paste0(names(fam.prev),
                                                           "(", fam.prev, "%)"), collapse = "\n"), "\n")
  cat("\nand the most prevalent genera are:\n", paste(paste0(names(gen.prev),
                                                             "(", gen.prev, "%)"), collapse = "\n"), "\n\n")
  phy.abund <- sort(phy.abund, decr = T)
  phy.abund <- round(phy.abund[phy.abund >= 0.05] * 100, 2)
  fam.abund <- sort(fam.abund, decr = T)
  fam.abund <- round(fam.abund[fam.abund >= 0.05] * 100, 2)
  gen.abund <- sort(gen.abund, decr = T)
  gen.abund <- round(gen.abund[gen.abund >= 0.05] * 100, 2)
  cat("\nThe most abundant phyla are:\n", paste(paste0(names(phy.abund),
                                                       "(", phy.abund, "%)"), collapse = "\n"), "\n")
  cat("\nThe most abundant families are:\n", paste(paste0(names(fam.abund),
                                                          "(", fam.abund, "%)"), collapse = "\n"), "\n")
  cat("\nand the most abundant genera are:\n", paste(paste0(names(gen.abund),
                                                            "(", gen.abund, "%)"), collapse = "\n"), "\n\n")
  sink()
}





OLD.perform_sequence_stat_analysis <- function(data.obj, ann = "") {
  sink(paste0("Sequence_Analysis_Statistics_", ann, ".txt"))
  otu.tab <- data.obj$otu.tab
  # Sequencing depth
  otu.abund <- rowSums(otu.tab)
  sam.abund <- colSums(otu.tab)
  otu.prev <- rowSums(otu.tab != 0)/ncol(otu.tab)
  otu.abund <- otu.abund[otu.abund >= 1]
  sam.abund <- sam.abund[sam.abund >= 1]
  cat("This data set contains ", length(sam.abund), " samples after quality controls.\n")
  cat("16S rDNA targeted sequencing yields ", median(sam.abund), "reads/sample on average (range:",
      min(sam.abund), "-", max(sam.abund), ").\n")
  cat("Clustering of these 16S sequence tags produces ", sum(otu.abund >
                                                               0), " OTUs at 97% similarity level.\n")
  pdf(paste0("Sequence_Analysis_Statistics_", ann, ".pdf"), height = 5,
      width = 5)
  obj <- ggplot2::ggplot(data = data.frame(x = otu.abund), aes(x = x)) +
    geom_histogram(col = "black", fill = "gray") + ylab("Frequency") +
    xlab("Abundance(Total counts)") + scale_x_log10(breaks = c(1, 10,
                                                               100, 1000, 10000, 1e+05, 1e+05))
  print(obj)
  obj <- ggplot2::ggplot(data = data.frame(x = sam.abund), aes(x = x)) +
    geom_histogram(col = "black", fill = "gray") + ylab("Frequency") +
    xlab("Sequencing depth")
  print(obj)
  obj <- ggplot2::ggplot(data = data.frame(x = otu.prev), aes(x = x)) +
    ylab("Frequency") + xlab("Prevalence(Occurence frequency)") + geom_histogram(col = "black",
                                                                                 fill = "gray")
  print(obj)
  dev.off()
  phy.abund <- data.obj$abund.list[["Phylum"]]
  fam.abund <- data.obj$abund.list[["Family"]]
  gen.abund <- data.obj$abund.list[["Genus"]]
  phy.prev <- rowSums(phy.abund != 0)/ncol(phy.abund)
  fam.prev <- rowSums(fam.abund != 0)/ncol(phy.abund)
  gen.prev <- rowSums(gen.abund != 0)/ncol(phy.abund)
  phy.abund <- rowMeans(t(t(phy.abund)/sam.abund))
  fam.abund <- rowMeans(t(t(fam.abund)/sam.abund))
  gen.abund <- rowMeans(t(t(gen.abund)/sam.abund))
  cat("These OTUs belong to ", sum(phy.abund > 0), " phyla,", sum(fam.abund >
                                                                    0), " families and ", sum(gen.abund > 0), "genera.\n\n")
  phy.prev <- sort(phy.prev, decr = T)
  phy.prev <- round(phy.prev[phy.prev >= 0.9] * 100, 2)
  fam.prev <- sort(fam.prev, decr = T)
  fam.prev <- round(fam.prev[fam.prev >= 0.9] * 100, 2)
  gen.prev <- sort(gen.prev, decr = T)
  gen.prev <- round(gen.prev[gen.prev >= 0.9] * 100, 2)
  # Rev: 2017_02_19 ' ' -> '\n'
  cat("\nThe most prevalent phyla are:\n", paste(paste0(names(phy.prev),
                                                        "(", phy.prev, "%)"), collapse = "\n"), "\n")
  cat("\nThe most prevalent families are:\n", paste(paste0(names(fam.prev),
                                                           "(", fam.prev, "%)"), collapse = "\n"), "\n")
  cat("\nand the most prevalent genera are:\n", paste(paste0(names(gen.prev),
                                                             "(", gen.prev, "%)"), collapse = "\n"), "\n\n")
  phy.abund <- sort(phy.abund, decr = T)
  phy.abund <- round(phy.abund[phy.abund >= 0.05] * 100, 2)
  fam.abund <- sort(fam.abund, decr = T)
  fam.abund <- round(fam.abund[fam.abund >= 0.05] * 100, 2)
  gen.abund <- sort(gen.abund, decr = T)
  gen.abund <- round(gen.abund[gen.abund >= 0.05] * 100, 2)
  cat("\nThe most abundant phyla are:\n", paste(paste0(names(phy.abund),
                                                       "(", phy.abund, "%)"), collapse = "\n"), "\n")
  cat("\nThe most abundant families are:\n", paste(paste0(names(fam.abund),
                                                          "(", fam.abund, "%)"), collapse = "\n"), "\n")
  cat("\nand the most abundant genera are:\n", paste(paste0(names(gen.abund),
                                                            "(", gen.abund, "%)"), collapse = "\n"), "\n\n")
  sink()
}
