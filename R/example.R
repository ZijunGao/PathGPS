# metabolites example

#' latentMediatorExample
#'
#' @param X ...
#'
#' @export
latentMediatorExample = function(data = "metabolites", plotOnly = TRUE){
  if(data == "metabolites"){
    rm(list = ls())
    # hyper-parameters
    p = 50; pNull = 400; r = 5 # r = 5 for tSNE/UMAP + kmeans; r = 8 for spectral clustering
    parametersUV = list()
    parametersUV$rotate.method = "varimax"; parametersUV$rotate.lambda = 0; parametersUV$rotate.m = 1; parametersUV$n.subsample = 200; parametersUV$randomSeed = 318; parametersUV$similarity.truncate = TRUE; parametersUV$probUZero = 0.6; parametersUV$probVZero = 0.6

    parametersClustering = list();
    parametersClustering$visualization.method = "UMAP" # "tSNE", *"UMAP", NULL
    parametersClustering$clustering.method = "kmeans"  # *"kmeans", "spectral"
    parametersClustering$number.center = seq(6,8) #  *8, "label propagation", seq(6,8)
    parametersClustering$pValue.k = 100; parametersClustering$pValue.threshold = NULL; parametersClustering$rsnps.read = FALSE; parametersClustering$save.address = NULL; parametersClustering$randomSeed = 318; parametersClustering$genotypeWeight = 0.0; parametersClustering$betweenWeight = 0.5; parametersClustering$phenotypeWeight = 1; parametersClustering$data = "metabolites" # pValue.k = 150;

    # estimate U, V
    SR = latentMediator::columnEstimator(betas = "metabolites", p = p, pNull = pNull, r = r, V.subtract = TRUE, rotate = TRUE, U.project = TRUE, UV.truncate = TRUE, parameters = parametersUV, aggregate.method = "bootstrap") # "bootstrap"

    # construct graph
    g = latentMediator::EVToGraphFull(U = SR$UHat, V = SR$VHat, kU = 4, kV = 8, betweenWeight = parametersClustering$betweenWeight, phenotypeWeight = parametersClustering$phenotypeWeight, genotypeWeight = parametersClustering$genotypeWeight, weightMatrixBootstrap = SR$weightMatrixBootstrap) #; kU = 3, kV = 6; SR$weightMatrixBootstra6

    # tSNE clustering plot
    pl = latentMediator::clusterPlot(g, maximalSimilarity = 1, parameters = parametersClustering)
    if(parametersClustering$visualization.method %in% c("tSNE", "UMAP")){
      ggplot2::ggsave("/Users/zijungao/Desktop/test2.pdf", pl$pl, width = 12, height = 8) # to be deleted
    }
  } else if (data == "biobankSub"){
    rm(list = ls())
    # hyper-parameters
    p = 1200; pNull = 250; r = 8 # r = 8 for tSNE/UMAP + kmeans; 11 for spectral

    parametersUV = list()
    parametersUV$rotate.method = "varimax"; parametersUV$rotate.lambda = 0; parametersUV$rotate.m = 1; parametersUV$n.subsample = 100; parametersUV$randomSeed = 318;
    parametersUV$similarity.truncate = TRUE; parametersUV$probUZero = 0.9; parametersUV$probVZero = 0.9

    parametersClustering = list();
    parametersClustering$visualization.method = "UMAP" # "tSNE", *"UMAP", NULL
    parametersClustering$tSNE.perplexity = 6 # *3
    parametersClustering$UMAP.n_neighbors = 8 # *5
    parametersClustering$clustering.method = "kmeans"  # *"kmeans", "spectral"
    parametersClustering$number.center = seq(9,12) # seq(9,12) for tSNE/UMAP + kmeans; seq(9,11) for spectral clustering
    parametersClustering$pValue.k = 100; parametersClustering$pValue.threshold = 0.05/888; parametersClustering$rsnps.read = FALSE; parametersClustering$save.address = NULL; parametersClustering$randomSeed = 318; parametersClustering$genotypeWeight = 0.05; parametersClustering$betweenWeight = 0.5; parametersClustering$phenotypeWeight = 1; parametersClustering$data = "biobankSub" # pValue.k = 150; weights: (0.05,0.5,1)

    # estimate U, V
    startTime = proc.time()
    SR = latentMediator::columnEstimator(betas = "biobankSub", p = p, pNull = pNull, pSample = floor(1.5 * p), r = r, V.subtract = TRUE, rotate = TRUE, U.project = TRUE, UV.truncate = TRUE, parameters = parametersUV, aggregate.method = "bootstrap") # "bootstrap"
    endTime = proc.time(); print("SR"); print(endTime[3] - startTime[3])

    # saveRDS(SR, file = "/Users/zijungao/Desktop/biobankTemp/UMAPSR.rds") # UMAPSR.rds/spectralSR.rds
    # SR = readRDS(file = "/Users/zijungao/Desktop/biobankTemp/UMAPSR.rds")

    # tSNE clustering plot
    # pl = latentMediator::clusterPlot(g, maximalSimilarity = 1, tSNE.perplexity = 5, clusteringMethod = "kmeans", clustering.centers = "label propagation", pValue.k = 150, pValue.threshold = NULL, rsnps.read = FALSE, saveAddress = NULL, randomSeed = 100)

    # construct graph
    startTime = proc.time()
    g = latentMediator::EVToGraphFull(U = SR$UHat, V = SR$VHat, kU = 3, kV = 6, betweenWeight = parametersClustering$betweenWeight, phenotypeWeight = parametersClustering$phenotypeWeight, genotypeWeight = parametersClustering$genotypeWeight, weightMatrixBootstrap = SR$weightMatrixBootstrap) # kU = 3, kV = 6/12; SR$weightMatrixBootstrap
    endTime = proc.time(); print("g"); print(endTime[3] - startTime[3])

    # tSNE clustering plot
    startTime = proc.time()
    parametersClustering$randomSeed = 1 # 1 for tSNE/UMAP + kmeans; 100 for spectral clustering
    pl = latentMediator::clusterPlot(g, maximalSimilarity = 1, parameters = parametersClustering)
    if(parametersClustering$visualization.method %in% c("tSNE", "UMAP")){
      ggplot2::ggsave("/Users/zijungao/Desktop/test3.pdf", pl$pl, width = 12, height = 8) # to be deleted
    }
    endTime = proc.time(); print("pl"); print(endTime[3] - startTime[3])
  } else if (data == "biobank"){
    rm(list = ls())
    # hyper-parameters
    p = 1200; pNull = 250; r = 20 # p = 1200; pNull = 250; r = 12 for tSNE/UMAP + kmeans; r = 20 for spectral clustering

    parametersUV = list()
    parametersUV$rotate.method = "varimax"; parametersUV$rotate.lambda = 0; parametersUV$rotate.m = 1; parametersUV$n.subsample = 100; parametersUV$randomSeed = 318;
    parametersUV$similarity.truncate = TRUE; parametersUV$probUZero = 0.9; parametersUV$probVZero = 0.9

    parametersClustering = list();
    parametersClustering$visualization.method = "tSNE" # "tSNE", *"UMAP", NULL
    parametersClustering$tSNE.perplexity = 6 # *3
    parametersClustering$UMAP.n_neighbors = 8 # *5
    parametersClustering$clustering.method = "kmeans"  # *"kmeans", "spectral"
    parametersClustering$number.center = seq(16,24, by = 2) # seq(9,12) for tSNE/UMAP + kmeans; seq(9,11) for spectral clustering
    parametersClustering$pValue.k = 100; parametersClustering$pValue.threshold = 0.05/2717; parametersClustering$rsnps.read = FALSE; parametersClustering$save.address = NULL; parametersClustering$randomSeed = 318; parametersClustering$genotypeWeight = 0.05; parametersClustering$betweenWeight = 0.5; parametersClustering$phenotypeWeight = 1; parametersClustering$data = "biobank" # pValue.k = 150; weights: (0.05,0.5,1)

    # estimate U, V
    # similarity.truncate = FALSE ~ 2min; similarity.truncate = TRUE ~ 6min
    startTime = proc.time()
    SR = latentMediator::columnEstimator(betas = "biobank", p = p, pNull = pNull, pSample = floor(1.5 * p), r = r, V.subtract = TRUE, rotate = TRUE, U.project = TRUE, UV.truncate = TRUE, parameters = parametersUV, aggregate.method = "bootstrap") # "bootstrap"
    endTime = proc.time(); print("SR"); print(endTime[3] - startTime[3])

    # saveRDS(SR, file = "/Users/zijungao/Desktop/biobankTemp/spectralSRFull.rds") # UMAPSRFull.rds/spectralSRFull.rds
    SR = readRDS(file = "/Users/zijungao/Desktop/biobankTemp/UMAPSRFull.rds")

    # tSNE clustering plot
    # pl = latentMediator::clusterPlot(g, maximalSimilarity = 1, tSNE.perplexity = 5, clusteringMethod = "kmeans", clustering.centers = "label propagation", pValue.k = 150, pValue.threshold = NULL, rsnps.read = FALSE, saveAddress = NULL, randomSeed = 100)

    # construct graph
    startTime = proc.time()
    g = latentMediator::EVToGraphFull(U = SR$UHat, V = SR$VHat, kU = 3, kV = 12, betweenWeight = parametersClustering$betweenWeight, phenotypeWeight = parametersClustering$phenotypeWeight, genotypeWeight = parametersClustering$genotypeWeight, weightMatrixBootstrap = SR$weightMatrixBootstrap) # kU = 3, kV = 6/12; SR$weightMatrixBootstrap
    endTime = proc.time(); print("g"); print(endTime[3] - startTime[3])

    # tSNE clustering plot
    startTime = proc.time()
    parametersClustering$randomSeed = 1 # 1 for tSNE/UMAP + kmeans; 100 for spectral clustering
    parametersClustering$tSNE.perplexity = 3
    parametersClustering$UMAP.n_neighbors = 15
    pl = latentMediator::clusterPlot(g, maximalSimilarity = 1, parameters = parametersClustering)
    if(parametersClustering$clustering.method %in% c("kmeans")){
      ggplot2::ggsave("/Users/zijungao/Desktop/test.pdf", pl$pl, width = 12, height = 8) # to be deleted
    }
    write.table(pl$result, "/Users/zijungao/Desktop/test.tsv", sep = "\t")
    endTime = proc.time(); print("pl"); print(endTime[3] - startTime[3])

  }
  result = list()
  if(!plotOnly){result$summary = pl$result; result$pl = pl$pl; result$VHat = SR$VHat; result$UHat = SR$UHat}
  return(result)
}

