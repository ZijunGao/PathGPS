# metabolites example

#' latentMediatorExample
#'
#' @param X ...
#'
#' @export
latentMediatorExample = function(data = "metabolites", plotOnly = TRUE, clustering.method = "UMAP"){
  if(data == "metabolites"){
    # rm(list = ls())
    suppressWarnings(RNGkind(sample.kind = "Rounding"))
    # hyper-parameters
    p = 50; p0 = 400; r = 6
    parameters = list()
    parameters$sprs.transform.m = 1; parameters$sprs.transform.lambda = 0
    parameters$U.prop.zero = 0.6; parameters$V.prop.zero = 0.6
    parameters$aggregate.method = "bootstrap"; parameters$n.aggregate = 300
    parameters$randomSeed = 100

    parameters$kU = 3; parameters$kV = 20
    parameters$gg.weight=0.1; parameters$gp.weight=1; parameters$pp.weight=1
    parameters$clustering.method = clustering.method
    parameters$max.coappearance = 1
    parameters$number.center = seq(6,7)
    parameters$tSNE.perplexity = NULL
    parameters$UMAP.n_neighbors = NULL
    parameters$clustering.randomSeed = 318
    parameters$coappearance.truncate  = 0.05

    parameters$pValue.k = 100
    parameters$pValue.threshold = NULL
    parameters$data = "metabolites"

    # estimation
    SR = pathGPS(beta = "metabolites", p = p, p0 = p0, r = r, envr.subtract = TRUE, sprs.transform = TRUE, aggregate = TRUE, parameters = parameters, returnFull = TRUE)
    # plot
    pl = clusterPlot(model = SR, pValue = "metabolites", names = TRUE)
    print(pl)
  # if(parameters$clustering.method %in% c("tSNE", "UMAP")){
  #     ggplot2::ggsave("/Users/zijungao/Desktop/test2.pdf", pl, width = 8, height = 8) # to be deleted
  #   }
  #   write.csv(SR$cluster, "/Users/zijungao/Desktop/test2.csv", row.names = FALSE)
    return(SR)
  }
  if (data == "biobank"){
    # rm(list = ls()); suppressWarnings(RNGkind(sample.kind = "Rounding"))
    # hyper-parameters
    # p = 1200; p0 = 250; r = 12 # p = 1200; p0 = 250; r = 12 for tSNE/UMAP + kmeans; r = 20 for spectral clustering

    # parametersUV = list()
    # parametersUV$rotate.method = "varimax"; parametersUV$rotate.lambda = 0; parametersUV$rotate.m = 1; parametersUV$n.aggregate = 100; parametersUV$randomSeed = 318; # parametersUV$n.subsample = 100
    # parametersUV$similarity.truncate = FALSE; parametersUV$U.prop.zero = 0.9; parametersUV$V.prop.zero = 0.9; parametersUV$aggregate.method = "bootstrap"
    #
    # parametersClustering = list();
    # parametersClustering$visualization.method = "UMAP" # "tSNE", *"UMAP", NULL
    # parametersClustering$tSNE.perplexity = 6 # *3
    # parametersClustering$UMAP.n_neighbors = 8 # *5
    # parametersClustering$clustering.method = "kmeans"  # *"kmeans", "spectral"
    # parametersClustering$number.center = seq(12,15, by = 1) # seq(9,12) for tSNE/UMAP + kmeans; seq(9,11) for spectral clustering; seq(6,24, by = 2)
    # parametersClustering$pValue.k = 100; parametersClustering$pValue.threshold = 0.05/2717; parametersClustering$rsnps.read = FALSE; parametersClustering$save.address = NULL; parametersClustering$randomSeed = 318; parametersClustering$gg.weight = 0.05; parametersClustering$gp.weight = 1; parametersClustering$pp.weight = 1; parametersClustering$data = "biobank" # pValue.k = 150; weights: (0.05,?,1)

    p = 1200; p0 = 250; r = 12 # p = 1200; p0 = 250; r = 12 for tSNE/UMAP + kmeans; r = 20 for spectral clustering

    parameters = list()
    parameters$sprs.transform.m = 1; parameters$sprs.transform.lambda = 0
    parameters$U.prop.zero = 0.9; parameters$V.prop.zero = 0.9
    parameters$aggregate.method = "bootstrap"; parameters$n.aggregate = 200
    parameters$randomSeed = 318
    parameters$coappearance.truncate = 0.01

    parameters$kU = 3; parameters$kV = 12 # 10, 20
    parameters$gg.weight=0.05; parameters$gp.weight=1; parameters$pp.weight=1
    parameters$clustering.method = "UMAP"
    parameters$max.coappearance = 1
    parameters$number.center = seq(8,15, by = 1)
    parameters$tSNE.perplexity = 3
    parameters$UMAP.n_neighbors = 15
    parameters$clustering.randomSeed = 100

    parameters$pValue.k = 100
    parameters$pValue.threshold = 0.05/2717
    parameters$data = "biobank"

    # estimate U, V
    # similarity.truncate = FALSE ~ 2min; similarity.truncate = TRUE ~ 6min
    startTime = proc.time()
    # SR = latentMediator::pathGPS(beta = "biobank", p = p, p0 = p0, pSample = floor(1.5 * p), r = r, V.subtract = TRUE, rotate = TRUE, U.project = TRUE, UV.truncate = TRUE, parameters = parametersUV, aggregate.method = "bootstrap") # "bootstrap"; V.subtract = TRUE,
    # model.biobank = pathGPS(beta = "biobank", p = p, p0 = p0, r = r, envr.subtract = TRUE, sprs.transform = TRUE, aggregate = TRUE, parameters = parameters, returnFull = TRUE)
    model.biobank = pathGPS(beta = "biobank", UV = model.biobank0[c(5,6)], p = p, p0 = p0, r = r, envr.subtract = TRUE, sprs.transform = TRUE, aggregate = TRUE, parameters = parameters, returnFull = TRUE) # UV = model.biobank0[c(5,6)]
    endTime = proc.time(); print("SR"); print(endTime[3] - startTime[3])
    # saveRDS(SR, file = "/Users/zijungao/Desktop/biobankTemp/UMAPSR0720.rds") # UMAPSRFull.rds/spectralSRFull.rds/UMAPSRFull1500.rds/spectralSRFull1500.rds spectralRotateFull.rds  spectralSR0720.rds/UMAPSR0720.rds
    # SR = readRDS(file = "/Users/zijungao/Desktop/biobankTemp/UMAPSR0720.rds")

    # construct graph
    # startTime = proc.time()
    # g = latentMediator::EVToGraphFull(U = SR$UHat, V = SR$VHat, kU = 10, kV = 20, gp.weight = parametersClustering$gp.weight, pp.weight = parametersClustering$pp.weight, gg.weight = parametersClustering$gg.weight, weightMatrixBootstrap = NULL)
    # weightMatrix = as.matrix(igraph::get.adjacency(g, type=c("both"), attr= "weight", names=TRUE)); diag(weightMatrix) = 0
    # weightMatrixZero = weightMatrix; weightMatrixZero[weightMatrixZero <= 0.02] = 0 # 0.01 for UMAP
    # isolatedNode = colnames(weightMatrix)[which(apply(weightMatrixZero,2,sum) == 0)]
    # if(length(isolatedNode) > 0) {g = igraph::delete.vertices(g, isolatedNode)}

    # endTime = proc.time(); print("g"); print(endTime[3] - startTime[3])

    # tSNE clustering plot
    # startTime = proc.time()
    # parametersClustering$randomSeed = 100 # 123 for tSNE; 100 for UMAP; 100 for spectral clustering
    # parametersClustering$randomSeed2 = 3 # 5 for tSNE; 3 for UMAP; 4 for spectral clustering
    # parametersClustering$tSNE.perplexity = 3
    # parametersClustering$UMAP.n_neighbors = 7
    # parametersClustering$clustering.method = "kmeans"; parametersClustering$visualization.method = "UMAP"
    # pl = latentMediator::clusterPlot(g, maximalSimilarity = 1, parameters = parametersClustering, names = FALSE)
    # endTime = proc.time(); print("pl"); print(endTime[3] - startTime[3])

    pl.biobank = clusterPlot(model = model.biobank, pValue = "biobank", names = FALSE)
    if(parameters$clustering.method %in% c("UMAP","tSNE")){
      ggplot2::ggsave("/Users/zijungao/Desktop/test2.pdf", pl.biobank, width = 6, height = 9) # to be deleted; with names: width = 12, height = 8; without names: width = 9, height = 6
    }
    write.table(model.biobank$cluster, "/Users/zijungao/Desktop/test2.tsv", sep = "\t", row.names = FALSE) # pl$result[,-c("eQTL")]
    # write.csv(model.biobank$cluster, "/Users/zijungao/Desktop/test2.csv", row.names = FALSE) # pl$result[,-c("eQTL")]
    endTime = proc.time(); print("pl"); print(endTime[3] - startTime[3])

  }
  result = list()
  if(!plotOnly){result$summary = pl$result; result$pl = pl$pl; result$VHat = SR$VHat; result$UHat = SR$UHat}
  # saveRDS(pl, "/Users/zijungao/Dropbox/latent mediator/plot/source code/source data/tSNEBiobank.rds")
  return(result)
}


# only 3 clusters, one very large!
# clp = igraph::cluster_label_prop(g)
# too slow
# betweenness = igraph::cluster_edge_betweenness(g)
# 7 clusters: too few clusters
# fastGreedy = igraph::cluster_fast_greedy(g)
# 6 clusters
# leadingEigen = igraph::cluster_leading_eigen(g)
# 7 clusters
# louvain = igraph::cluster_louvain(g)
# GLPK is not available
# optimal = igraph::cluster_optimal(g)
# slow; 21 clusters, but ~10 has less than 5 nodes. Essentially 8 clusters
# spinglass = igraph::cluster_spinglass(g)
# 10 clusters; but a giant cluster (>50%)
# randomWalk = igraph::cluster_walktrap(g)

# playground
# temp = read.table(file = "./data/biobankFullBeta.txt", header = TRUE)
# description = read.table(file = "./data/biobankDescription3.txt", header = TRUE)
# temp = temp[,description$name]
# write.table(temp, file = "./data/biobankBeta.txt")
# test = read.table(file = "./data/biobankBeta.txt", header = TRUE)

# SR$pl = pl; saveRDS(SR, file = "/Users/zijungao/Dropbox/PathGPS/plot/source code/source data/metabolitesSpectral.rds")

