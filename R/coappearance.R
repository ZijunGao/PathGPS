# Extract the list of edges ordered according to co-appearance weights
#'
#' order each type of edges decreasingly in edge weights.
#' @param g a weighted graph object from igraph. By default, there are three types of edges: gene-gene, gene-trait, and trait-trait. Genes are recorded as "rsXXX", and traits do not start with "rs".
#' @param top how many edges for display.
#' @import igraph
#'
#' @examples
#' SR0 = readRDS(file = "/Users/zijungao/Dropbox/PathGPS/plot/source code/source data/metabolitesUMAP.rds")
#' coappearance.list(SR0$graph)
#'
#' @export
coappearance.list = function(g, top = 10){
  edge.weight = as.data.frame(cbind(E(g)$weight, get.edgelist(g)))
  colnames(edge.weight) = c("weight", "vertex1", "vertex2")
  edge.weight$weight = as.numeric(edge.weight$weight)
  edge.weight$type = factor((substring(edge.weight$vertex1, first = 1, last = 2) == "rs") + (substring(edge.weight$vertex2, first = 1, last = 2) == "rs"))
  levels(edge.weight$type) = c("trait-trait", "gene-trait", "gene-gene")

  result = list()
  # order each type of edges decreasingly in weights
  edge.weight.pp = edge.weight[edge.weight$type == "trait-trait",]
  result$trait.trait = edge.weight.pp[order(edge.weight.pp$weight, decreasing = T)[1:min(top, dim(edge.weight.pp)[1])],]

  edge.weight.gp = edge.weight[edge.weight$type == "gene-trait",]
  result$gene.trait = edge.weight.gp[order(edge.weight.gp$weight, decreasing = T)[1:min(top, dim(edge.weight.gp)[1])],]

  edge.weight.gg = edge.weight[edge.weight$type == "gene-gene",]
  result$gene.gene = edge.weight.gg[order(edge.weight.gg$weight, decreasing = T)[1:min(top, dim(edge.weight.gg)[1])],]

  return(result)
}

# index = which(get.edgelist(SR.biobank$graph)[,2] == "X4105_irnt")
# cbind(get.edgelist(SR.biobank$graph), E(SR.biobank$graph)$weight)[index,]


