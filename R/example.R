# metabolites example

#' latentMediatorExample
#'
#' @param X ...
#'
#' @export
latentMediatorExample = function(data = "metabolites", plotOnly = TRUE, clustering.method = "UMAP"){
  if(data == "metabolites"){
    # rm(list = ls())
    # hyper-parameters
    p = 50; p0 = 400; r = 4
    parameters = list()
    parameters$sprs.transform.m = 1; parameters$sprs.transform.lambda = 0
    parameters$U.prop.zero = 0.6; parameters$V.prop.zero = 0.6
    parameters$aggregate.method = "bootstrap"; parameters$n.aggregate = 300
    parameters$randomSeed = 100

    parameters$kU = 3; parameters$kV = 6
    parameters$gg.weight=0.1; parameters$gp.weight=0.5; parameters$pp.weight=1
    parameters$clustering.method = clustering.method
    parameters$max.coappearance = 1
    parameters$number.center = seq(6,7)
    parameters$tSNE.perplexity = NULL
    parameters$UMAP.n_neighbors = NULL
    parameters$clustering.randomSeed = 318
    parameters$coappearance.truncate  = NULL

    parameters$pValue.k = 150
    parameters$pValue.threshold = NULL
    parameters$data = "metabolites"

    # estimation
    SR = pathGPS(beta = "metabolites", UV = SR0[c(5,6)], p = p, p0 = p0, r = r, envr.subtract = TRUE, sprs.transform = TRUE, aggregate = TRUE, parameters = parameters, returnFull = TRUE)
    # plot
    pl = clusterPlot(model = SR, pValue = "metabolites", names = F)
    pl.name = clusterPlot(model = SR, pValue = "metabolites", names = T)
    # metabolites.coappearance = read.table("./data/metabolitesCoappearance.txt")
    # pl.metabolites = clusterPlot(model = SR, pValue = as.matrix(metabolites.coappearance), names = F)
    # pl.metabolites.name = clusterPlot(model = SR, pValue = as.matrix(metabolites.coappearance), names = T)
    # print(pl)
  if(parameters$clustering.method %in% c("tSNE", "UMAP")){
      ggplot2::ggsave("/Users/zijungao/Desktop/test2.pdf", pl.name, width = 8, height = 8) # to be deleted
    }
    write.csv(SR$cluster, "/Users/zijungao/Desktop/test2.csv", row.names = FALSE)
    # SR0 = readRDS(file = "/Users/zijungao/Dropbox/PathGPS/plot/source code/source data/metabolitesUMAP.rds")
    return(SR)
  }
  if (data == "biobank"){
    # rm(list = ls());
    # hyper-parameters

    p = 1200; p0 = 250; r = 12 # p = 1200; p0 = 250; r = 12 for tSNE/UMAP + kmeans; r = 20 for spectral clustering

    parameters = list()
    parameters$sprs.transform.m = 1; parameters$sprs.transform.lambda = 0
    parameters$U.prop.zero = 0.9; parameters$V.prop.zero = 0.9
    parameters$aggregate.method = "bootstrap"; parameters$n.aggregate = 200
    parameters$randomSeed = 318
    parameters$coappearance.truncate = 0.01

    parameters$kU = 3; parameters$kV = 12 # 12, 12, 20
    parameters$gg.weight=0.05; parameters$gp.weight=1; parameters$pp.weight=1
    parameters$clustering.method = "UMAP"
    parameters$max.coappearance = 1
    parameters$number.center = seq(8,15, by = 1)
    parameters$tSNE.perplexity = 6
    parameters$UMAP.n_neighbors = 15
    parameters$clustering.randomSeed = 100 # 100; 100; 123

    parameters$pValue.k = 350
    parameters$pValue.threshold = 0.05/2717
    parameters$data = "biobank"

    # estimate U, V
    model.biobank = pathGPS(beta = "biobank", UV = model.biobank0[c(5,6)], p = p, p0 = p0, r = r, envr.subtract = TRUE, sprs.transform = TRUE, aggregate = TRUE, parameters = parameters, returnFull = TRUE) # UV = NULL or record (model.biobank0[c(5,6)])

    pl.biobank = clusterPlot(model = model.biobank, pValue = "biobank", names = F)
    # biobank.coappearance = read.table("./data/biobankCoappearance.txt")
    # pl.biobank = clusterPlot(model = model.biobank, pValue = as.matrix(biobank.coappearance), names = F)
    # pl.biobank.name = clusterPlot(model = model.biobank, pValue = as.matrix(biobank.coappearance), names = T)
    if(parameters$clustering.method %in% c("UMAP","tSNE")){
      ggplot2::ggsave("/Users/zijungao/Desktop/test2.pdf", pl.biobank, width = 6, height = 6) # to be deleted; with names: width = 12, height = 8; UMAP: width = 6, height = 15; tSNE: width = 6, height = 6
    }
    # write.table(model.biobank$cluster, "/Users/zijungao/Desktop/test2.tsv", sep = "\t", row.names = FALSE) # pl$result[,-c("eQTL")]
    endTime=proc.time(); print("pl"); print(endTime[3] - startTime[3])
  }
  result = list()
  if(!plotOnly){result$summary = pl$result; result$pl = pl$pl; result$VHat = SR$VHat; result$UHat = SR$UHat}
  # model.biobank$pl = pl.biobank2
  # model.biobank$pl.name = pl.biobank
  # saveRDS(model.biobank, "/Users/zijungao/Dropbox/PathGPS/plot/source code/source data/biobankUMAP.rds")
  # model.biobank/model.biobank0 = readRDS("/Users/zijungao/Dropbox/PathGPS/plot/source code/source data/biobankUMAP.rds") # model.biobank recovers the result Sept. 15, 2021, the most up-to-date
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

# prev

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

# Mar. 22
# # metabolites
# create the co-appearance matrix
# metabolites.pvalue = data.matrix(matrix(0, length(SR$cluster$name), length(SR$cluster$name)))
# rownames(metabolites.pvalue)  = colnames(metabolites.pvalue) = SR$cluster$name
# for(i in 1:3){
#   metabolites.coappearance = coappearance.list(SR$graph, top = 10000)[[i]] # 10000 > max{#genes, #traits}
#   metabolites.pvalue[cbind(metabolites.coappearance$vertex1,metabolites.coappearance$vertex2)] = metabolites.coappearance$weight
# }
# metabolites.pvalue.transform = exp(-100*metabolites.pvalue) # hyperparameter: 100
# # UK Biobank
# create the coappearance matrix
# biobank.pvalue = data.matrix(matrix(0, length(model.biobank$cluster$name), length(model.biobank$cluster$name)))
# rownames(biobank.pvalue)  = colnames(biobank.pvalue) = model.biobank$cluster$name
# for(i in 1:3){
#   biobank.coappearance = coappearance.list(model.biobank$graph, top = 10000)[[i]] # 10000 > max{#genes, #traits}
#   biobank.pvalue[cbind(biobank.coappearance$vertex1,biobank.coappearance$vertex2)] = biobank.coappearance$weight
# }
# biobank.pvalue.transform = exp(-100*biobank.pvalue) # hyperparameter: 100


