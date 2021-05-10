# metabolites example

#' latentMediatorExample
#'
#' @param X ...
#'
#' @export
latentMediatorExample = function(data = "metabolites"){
  if(data == "metabolites"){
    # hyper-parameters
    p = 50; p2 = 400; r = 5; m = 100
    # estimate U, V
    SR = latentMediator::columnEstimator(betas = "metabolites", p = p, pNull = p2, pSample = p*2, r = r, m = m, V.subtract = TRUE, V.rotate = TRUE, V.rotate.m = 1, U.project = TRUE, randomSeed = 318)
    # construct graph
    g = latentMediator::EVToGraphFull(U = SR$UHat, V = SR$VHat, kU = 3, kV = 6, betweenWeight = 0.5, phenotypeWeight = 1, genotypeWeight = 0)
    # tSNE clustering plot
    pl = latentMediator::tSNEPlot(g, maximalSimilarity = 1, tSNE.perplexity = 5, clusteringMethod = "kmeans", clustering.centers = "label propagation", pValue.k = 150, pValue.threshold = NULL, rsnps.read = FALSE, saveAddress = NULL, randomSeed = 100)
  }
  result = list(); result$pl = pl; result$VHat = SR$VHat; result$UHat = SR$UHat;
  return(result)
}
