# visualization

#' tSNEPlot
#'
#' tSNE plot and kmeans clusters
#'
#' @param X covariate matrix to predict conditional densities at; each row represents an observation vector.
#'
#' @export
tSNEPlot = function(g, maximalSimilarity = 1, tSNE.perplexity = 5, clusteringMethod = "kmeans", clustering.centers = "label propagation", pValue.k = 150, pValue.threshold = NULL, rsnps.read = FALSE, saveAddress = NULL, randomSeed = 100){
  # tSNE
  weightMatrix = as.matrix(igraph::get.adjacency(g, type=c("both"), attr= "weight", names=TRUE)); diag(weightMatrix) = maximalSimilarity
  distanceMatrix = maximalSimilarity - weightMatrix
  if(is.null(colnames(distanceMatrix))){stop("please input names of genotypes and phenotypes")}
  set.seed(randomSeed)
  tsne = Rtsne::Rtsne(X = distanceMatrix, is_distance = TRUE, perplexity = tSNE.perplexity)

  # prepare data
  if(clusteringMethod == "kmeans"){
    if(clustering.centers == "label propagation"){
      set.seed(randomSeed)
      clp = igraph::cluster_label_prop(g)
      clustering.centers = length(unique(igraph::membership(clp)))
    }
    set.seed(randomSeed)
    clusteringMembership = kmeans(tsne$Y, centers = clustering.centers)$cluster
  } else if(clusteringMethod == "label propagation"){
    set.seed(randomSeed)
    clp = igraph::cluster_label_prop(g)
    clusteringMembership = igraph::membership(clp)
  }
  data <- data.frame(tsne.x = tsne$Y[, 1],
                     tsne.y = tsne$Y[, 2],
                     group = as.factor(clusteringMembership),
                     name = colnames(distanceMatrix))
  # data = data.frame(tsne.x = tsne$Y[, 1], tsne.y = tsne$Y[, 2], group = as.factor(clusteringMembership), name = colnames(distanceMatrix))
  data = data.table::data.table(data)
  data[, genotype := grepl("^rs", name)]
  data[, size := 3 - genotype]

  # get neighbor genes for SNPs
  if(rsnps.read){
    httr::set_config(httr::config(http_version = 0))
    rsnps_results = rsnps::ncbi_snp_query(data[data$genotype, ]$name)
  } else {
    rsnps_results = read.table("./data/metabolitesRsnps.txt", header = TRUE)
  }
  eQTL_results = read.table("./data/metabolitesGeneEQTL.txt", header = TRUE) # TODO: replace as the source code
  rsnps_results = merge(rsnps_results, eQTL_results, by.x = "rsid", by.y = "rsid", all.x = TRUE)
  rsnps_results$gene[which((rsnps_results$gene != "") * (rsnps_results$geneEQTL != "") == 1)] = paste(rsnps_results$gene[which((rsnps_results$gene != "") * (rsnps_results$geneEQTL != "") == 1)], rsnps_results$geneEQTL[which((rsnps_results$gene != "") * (rsnps_results$geneEQTL != "") == 1)], sep = "/")
  rsnps_results$gene[which((rsnps_results$gene == "") * (rsnps_results$geneEQTL != "") == 1)] = paste(rsnps_results$gene[which((rsnps_results$gene == "") * (rsnps_results$geneEQTL != "") == 1)], rsnps_results$geneEQTL[which((rsnps_results$gene == "") * (rsnps_results$geneEQTL != "") == 1)], sep = "")
  rsnps_results$geneEQTL = NULL

  data = merge(data, rsnps_results[, c("rsid", "gene")], by.x = "name", by.y = "rsid", all.x = TRUE)
  data[, plot.name := name]
  data[gene != "", plot.name := gene]
  data[genotype == TRUE, plot.name := paste0("italic('", plot.name, "')")]

  # form edges between phenotypes and genotypes
  pValue = data.matrix(read.table("./data/metabolitesPValue.txt"))
  edge = reshape2::melt(pValue <= sort(pValue, decreasing = FALSE)[pValue.k])
  if(is.null(pValue.threshold)){pValue.threshold = 0.05/dim(pValue)[2]}
  edge2 = reshape2::melt(pValue <= pmin(pValue.threshold, apply(pValue, 1, min)))
  edge$value = (pmax(edge$value, edge2$value) >= 1)
  edge$pValue = reshape2::melt(pValue)$value
  edge = subset(edge, value)
  names(edge)[1:2] = c("genotype", "phenotype")
  edge = merge(edge, data[, c("name", "tsne.x", "tsne.y")], by.x = "genotype", by.y = "name")
  names(edge)[5:6] = c("genotype.x", "genotype.y")
  edge = merge(edge, data[, c("name", "tsne.x", "tsne.y")], by.x = "phenotype", by.y = "name")
  names(edge)[7:8] = c("phenotype.x", "phenotype.y")

  # plot
  pl = ggplot2::ggplot(data) + ggplot2::aes(x = tsne.x, y = tsne.y) + ggplot2::geom_point(ggplot2::aes(size = size, shape = genotype, color = group)) + ggrepel::geom_text_repel(parse = TRUE, alpha = 0.7, max.overlaps = 20, ggplot2::aes(size = size, label = plot.name))
  pl$layers = c(ggplot2::geom_segment(ggplot2::aes(x = genotype.x, y = genotype.y, xend = phenotype.x, yend = phenotype.y), size = 0.5-0.5*log(edge$pValue)/max(-log(edge$pValue)), data = edge, color = "gray", alpha = 0.3-0.5 * log(edge$pValue)/max(-log(edge$pValue))), pl$layers)
  pl = pl + ggplot2::theme_void() + ggplot2::scale_size_identity() + ggplot2::theme(legend.position = "none")
  if(!is.null(saveAddress)){
    ggplot2::ggsave("/Users/zijungao/Desktop/ggplot-visualizationSimple.pdf", pl, width = 12, height = 8) # ggplot-visualization.pdf/ggplot-visualizationVHatRotate.pdf
  } else {return(pl)}
}
