# visualization

#' clusterPlot
#'
#' clustering and visualization
#'
#' @param g ...
#' @param parameters: visualization.method, tSNE.perplexity, clustering.method, cluster.centers, pValue.k, pValue.threshold, rsnps.read, saveAddress, randomSeed, UMAP.n_neighbors, genotypeWeight, betweenWeight, phenotypeWeight
#'
#' @export
clusterPlot = function(g, maximalSimilarity = 1, parameters = list()){
  # preprocess
  if(is.null(parameters$randomSeed)){parameters$randomSeed = 100}
  if(is.null(parameters$visualization.method)){parameters$visualization.method = "tSNE"}
  if(is.null(parameters$clustering.method)){parameters$clustering.method = "kmeans"}
  if(parameters$clustering.method == "cmeans"){stop("cmeans clustering is currently not available")}

  weightMatrix = as.matrix(igraph::get.adjacency(g, type=c("both"), attr= "weight", names=TRUE)); diag(weightMatrix) = 0
  distanceMatrix = maximalSimilarity - weightMatrix; diag(distanceMatrix) = 0
  if(is.null(colnames(distanceMatrix))){stop("please input names of genotypes and phenotypes")}
  if(parameters$number.center[1] == "label propagation"){
    set.seed(parameters$randomSeed)
    clp = igraph::cluster_label_prop(g)
    parameters$number.center = length(unique(igraph::membership(clp)))
  }

  if(parameters$clustering.method == "spectral"){
    set.seed(parameters$randomSeed)
    kernMatrix = weightMatrix;
    if(parameters$data %in% c("metabolites", "biobankSub")){
      diag(kernMatrix) = 0.2
    } else if(parameters$data %in% c("biobank")){
      diag(kernMatrix) = 0.1
    }
    clusteringMembership = myClustering(kernMatrix = kernMatrix, weightMatrix = weightMatrix, method = "spectral", parameters = parameters)$membership
    D = apply(kernMatrix, 1, sum)
    spectral.embedding = eigen(diag(sqrt(1/D)) %*%  kernMatrix %*% diag(sqrt(1/D)))$vectors
    data <- data.frame(tsne.x = spectral.embedding[, 1],
                       tsne.y = spectral.embedding[, 2],
                         group = as.factor(clusteringMembership),
                         name = colnames(distanceMatrix))
  } else if (parameters$clustering.method == "hierarchical"){
    stop("hierarchical clustering can not be visualized")
    clusteringMembership = myClustering(distMatrix = distanceMatrix, weightMatrix = weightMatrix, method = "hierarchical", parameters = parameters)$membership
  } else if(parameters$visualization.method == "tSNE"){
    set.seed(parameters$randomSeed)
    if(is.null(parameters$tSNE.perplexity)){
      parameters$tSNE.perplexity = floor(min(max(5, sqrt(dim(distanceMatrix)[1])), (dim(distanceMatrix)[1]-1)/3))
    }
    tsne = Rtsne::Rtsne(X = distanceMatrix, is_distance = TRUE, perplexity = parameters$tSNE.perplexity)
    rownames(tsne$Y) = rownames(distanceMatrix)
    if(parameters$clustering.method == "label propagation"){
      set.seed(parameters$randomSeed)
      clp = igraph::cluster_label_prop(g)
      clusteringMembership = igraph::membership(clp)
    }
    if(parameters$clustering.method %in% c("kmeans", "cmeans")){
      clusteringMembership = myClustering(X = tsne$Y, weightMatrix = weightMatrix, method = parameters$clustering.method, parameters = parameters)$membership
    }
    data <- data.frame(tsne.x = tsne$Y[, 1],
                       tsne.y = tsne$Y[, 2],
                       group = as.factor(clusteringMembership),
                       name = colnames(distanceMatrix))
  } else if(parameters$visualization.method == "UMAP"){
    umap.parameters = umap::umap.defaults
    umap.parameters$n_components = 2; umap.parameters$input = "dist"
    if(!is.null(parameters$UMAP.n_neighbors)){
      umap.parameters$n_neighbors = parameters$UMAP.n_neighbors
    } else {umap.parameters$n_neighbors = floor(min(max(3, sqrt(dim(distanceMatrix)[1])), 15))}
    set.seed(parameters$randomSeed)
    umapEmbedding = umap::umap(d = distanceMatrix, config = umap.parameters, method = "naive")
    if(parameters$clustering.method == "label propagation"){
      set.seed(parameters$randomSeed)
      clp = igraph::cluster_label_prop(g)
      clusteringMembership = igraph::membership(clp)
    }
    if(parameters$clustering.method %in% c("kmeans", "cmeans")){
      clusteringMembership = myClustering(X =  umapEmbedding$layout, weightMatrix = weightMatrix, method = parameters$clustering.method, parameters = parameters)$membership
    }
    data <- data.frame(tsne.x = umapEmbedding$layout[, 1],
                       tsne.y = umapEmbedding$layout[, 2],
                       group = as.factor(clusteringMembership),
                       name = colnames(distanceMatrix))
  }

  # visualization
  # data <- data.frame(tsne.x = tsne$Y[, 1],
  #                    tsne.y = tsne$Y[, 2],
  #                    group = as.factor(clusteringMembership),
  #                    name = colnames(distanceMatrix))
  # data = data.frame(tsne.x = tsne$Y[, 1], tsne.y = tsne$Y[, 2], group = as.factor(clusteringMembership), name = colnames(distanceMatrix))
  data = data.table::data.table(data)
  data[, genotype := grepl("^rs", name)]
  data[, size := 3 - genotype]

  # get neighbor genes for SNPs
  if(parameters$rsnps.read){
    httr::set_config(httr::config(http_version = 0))
    rsnps_results = rsnps::ncbi_snp_query(data[data$genotype, ]$name)
  } else {
    if(parameters$data == "metabolites"){
      rsnps_results = read.table("./data/metabolitesRsnps.txt", header = TRUE)
    } else if(parameters$data == "biobankSub"){
      rsnps_results = read.table("./data/biobank888Rsnps.txt", header = TRUE)
    } else if(parameters$data == "biobank"){
      # TODO!!!
      rsnps_results = read.table("./data/biobank888Rsnps.txt", header = TRUE)
    }
  }
  if(parameters$data == "metabolites"){
    eQTL_results = read.table("./data/metabolitesGeneEQTL.txt", header = TRUE) # TODO: replace as the source code
  } else if (parameters$data == "biobankSub"){
    eQTL_results = read.table("./data/metabolitesGeneEQTL.txt", header = TRUE) # TODO: replace by biobank's version
  } else if (parameters$data == "biobank"){
    eQTL_results = read.table("./data/metabolitesGeneEQTL.txt", header = TRUE) # TODO: replace by biobank's version
  }
  if(parameters$data == "biobankSub"){
    description = read.table("./data/biobankSubDescription.txt", header = TRUE)
  } else if (parameters$data == "biobank"){
    description = read.table("./data/biobankDescription.txt", header = TRUE)
  }
  rsnps_results = merge(rsnps_results, eQTL_results, by.x = "rsid", by.y = "rsid", all.x = TRUE)
  rsnps_results$gene[which((rsnps_results$gene != "") * (rsnps_results$geneEQTL != "") == 1)] = paste(rsnps_results$gene[which((rsnps_results$gene != "") * (rsnps_results$geneEQTL != "") == 1)], rsnps_results$geneEQTL[which((rsnps_results$gene != "") * (rsnps_results$geneEQTL != "") == 1)], sep = "/")
  rsnps_results$gene[which((rsnps_results$gene == "") * (rsnps_results$geneEQTL != "") == 1)] = paste(rsnps_results$gene[which((rsnps_results$gene == "") * (rsnps_results$geneEQTL != "") == 1)], rsnps_results$geneEQTL[which((rsnps_results$gene == "") * (rsnps_results$geneEQTL != "") == 1)], sep = "")
  rsnps_results$geneEQTL = NULL

  if(parameters$clustering.method == "spectral"){
    data$position = seq(1, dim(kernMatrix)[1])
    data = merge(data, rsnps_results[, c("rsid", "gene")], by.x = "name", by.y = "rsid", all.x = TRUE)
    data[, plot.name := name]
    data[gene != "", plot.name := gene]
    colnames(kernMatrix)[data$position] = data$plot.name
    rownames(kernMatrix)[data$position] = data$plot.name
    clusteringSize = (table(clusteringMembership) + runif(table(clusteringMembership), 0, 0.01))[clusteringMembership]
    # genotype-phenotype cluster
    heatmap(kernMatrix[order(clusteringSize), order(clusteringSize)], Colv = NA, Rowv = NA, symm = TRUE, add.expr = abline(h=c(0.5,0.5+cumsum(sort(table(clusteringMembership), decreasing = FALSE))), v=c(0.5,0.5+cumsum(sort(table(clusteringMembership), decreasing = FALSE))), lty = 2, lwd = 1), ColSideColors = c("black", "red")[1+data$genotype[order(data$position)[order(clusteringSize)]]], RowSideColors = c("black", "red")[1+data$genotype[order(data$position)[order(clusteringSize)]]])
    # bi-cluster
    if(parameters$data == "metabolites"){
      pValue = data.matrix(read.table("./data/metabolitesPValue.txt"))
    } else if(parameters$data == "biobankSub"){
      pValue = data.matrix(read.table("./data/biobankSubPValue.txt"))
    } else if(parameters$data == "biobank"){
      pValue = data.matrix(read.table("./data/biobankPValue.txt"))
    }
    pValueBicluster = pValue[data$name[data$genotype], data$name[!data$genotype]]
    rownames(pValueBicluster) = data$plot.name[data$genotype]
    heatmap(-log(pValueBicluster)[order(clusteringSize[data$position[data$genotype]]), order(clusteringSize[data$position[!data$genotype]])], Colv = NA, Rowv = NA, symm = FALSE, add.expr = abline(h=c(0.5,0.5+cumsum(pmax(0.4,aggregate(data$genotype, by = list(data$group), sum)$x[order(table(clusteringMembership))]))), v=c(0.5,0.5+cumsum(pmax(0.4, aggregate(1-data$genotype, by = list(data$group), sum)$x[order(table(clusteringMembership))]))), lty = 2, lwd = 1)) # duplicated columns/rows
    # clustering
    clusteringMembership = rep(seq(1,length(table(clusteringMembership))), sort(table(clusteringMembership)))
    names(clusteringMembership) = rownames(kernMatrix)[order(clusteringSize)]

    result = data[,c("plot.name", "group", "genotype", "position")]
    colnames(result)[c(1,2)] = c("name", "cluster")
    result$cluster[order(result$position)[order(clusteringSize)]] = rep(seq(1,length(table(clusteringMembership))), sort(table(clusteringMembership)))
    description$name[substr(description$name,1,1) %in% as.character(seq(0,9))] = paste("X", description$name[substr(description$name,1,1) %in% as.character(seq(0,9))], sep = "") # ???
    result = merge(result, description, by.x = "name", by.y = "name", all.x = TRUE)
    result = result[order(result$cluster),]; result$position = NULL
    result2 = list(); result2$result = result; result2$pl = NA
    return(result2)
  }

  data = merge(data, rsnps_results[, c("rsid", "gene")], by.x = "name", by.y = "rsid", all.x = TRUE)
  data[, plot.name := name]
  data[gene != "", plot.name := gene]
  data[genotype == TRUE, plot.name := paste0("italic('", plot.name, "')")]
  # form edges between phenotypes and genotypes
  if(parameters$data == "metabolites"){
    pValue = data.matrix(read.table("./data/metabolitesPValue.txt"))
  } else if(parameters$data == "biobankSub"){
    pValue = data.matrix(read.table("./data/biobankSubPValue.txt"))
  } else if(parameters$data == "biobank"){
    pValue = data.matrix(read.table("./data/biobankPValue.txt"))
  }
  edge = reshape2::melt(pValue <= sort(pValue, decreasing = FALSE)[parameters$pValue.k])
  if(is.null(parameters$pValue.threshold)){parameters$pValue.threshold = 0.05/dim(pValue)[2]}
  edge2 = reshape2::melt(pValue <= pmin(parameters$pValue.threshold, apply(pValue, 1, min)))
  edge$value = (pmax(edge$value, edge2$value) >= 1)
  edge$pValue = reshape2::melt(pValue)$value
  edge = subset(edge, value)
  names(edge)[1:2] = c("genotype", "phenotype")
  edge = merge(edge, data[, c("name", "tsne.x", "tsne.y")], by.x = "genotype", by.y = "name")
  names(edge)[5:6] = c("genotype.x", "genotype.y")
  edge = merge(edge, data[, c("name", "tsne.x", "tsne.y")], by.x = "phenotype", by.y = "name")
  names(edge)[7:8] = c("phenotype.x", "phenotype.y")

  # plot
  if(parameters$data %in% c("metabolites", "biobankSub")){
    pl = ggplot2::ggplot(data) + ggplot2::aes(x = tsne.x, y = tsne.y) + ggplot2::geom_point(ggplot2::aes(size = size, shape = genotype, color = group)) + ggrepel::geom_text_repel(parse = TRUE, alpha = 0.7, max.overlaps = 20, ggplot2::aes(size = size, label = plot.name))
  } else if (parameters$data %in% c("biobank")){
    data$size = data$size - 0.5
    pl = ggplot2::ggplot(data) + ggplot2::aes(x = tsne.x, y = tsne.y) + ggplot2::geom_point(ggplot2::aes(size = size, shape = genotype, color = group)) + ggrepel::geom_text_repel(parse = TRUE, alpha = 0.7, max.overlaps = 30, ggplot2::aes(size = size, label = plot.name))
  }
  pl$layers = c(ggplot2::geom_segment(ggplot2::aes(x = genotype.x, y = genotype.y, xend = phenotype.x, yend = phenotype.y), size = 0.5-0.5*log(edge$pValue)/max(-log(edge$pValue)), data = edge, color = "gray", alpha = 0.3-0.5 * log(edge$pValue)/max(-log(edge$pValue))), pl$layers)
  pl = pl + ggplot2::theme_void() + ggplot2::scale_size_identity() + ggplot2::theme(legend.position = "none")

  if(!is.null(parameters$save.address)){
    ggplot2::ggsave("/Users/zijungao/Desktop/ggplot-visualizationSimple.pdf", pl, width = 12, height = 8) # ggplot-visualization.pdf/ggplot-visualizationVHatRotate.pdf
  }

  result = data[,c("plot.name", "group", "genotype")]
  colnames(result)[c(1,2)] = c("name", "cluster")
  result$name[result$genotype] = substr(result$name[result$genotype], start = 9, stop = nchar(result$name[result$genotype])-2)
  description$name[substr(description$name,1,1) %in% as.character(seq(0,9))] = paste("X", description$name[substr(description$name,1,1) %in% as.character(seq(0,9))], sep = "") # ???
  if(parameters$data %in% c("biobankSub", "biobank")){
    result = merge(result, description, by.x = "name", by.y = "name", all.x = TRUE)}
  result = result[order(result$cluster),]
  result2 = list(); result2$result = result; result2$pl = pl
  return(result2)
}
