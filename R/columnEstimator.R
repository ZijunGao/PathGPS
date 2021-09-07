# helper of metabolites

#' pathGPS
#'
#' estimate the row space of V (#mediators x #phenotypes) and the column space of U (#genotypes x #mediators)
#' @param X covariate matrix to predict conditional densities at; each row represents an observation vector.
#' @param beta full marginal association estimates, either a string or a matrix of dimension number of SNPs x number of phenotypes. If \code{beta} is "metabolites", the metabolites data is used; if \code{beta} is "biobank", the UK biobank data is used. If \code{beta} is a matrix, each row should stand for a SNP, each column should stand for a phenotype, and the value should be the marginal association estimate of the corresponding SNP, phenotype pair. SNPs should be ranked in decreasing order of their relevance to the set of phenotypes. The top \code{p} rows will be used as signal SNPs and the last \code{p0} rows will be used as noise SNPs.
#' @param p the number of signal SNPs.
#' @param p0 the number of noise SNPs.
#' @param r the number of eigen-vectors preserved in the truncated eigen-decomposition.
#' @param envr.subtract a logical value. If \code{TRUE}, environmental effects are estimated using noise SNPs and subtracted from genetic effects in signal SNPs. Default is \code{TRUE}.
#' @param sprs.transform a logical value. If \code{TRUE}, estimated eigen-vectors are linearly transformed to be sparser using \code{varimax} or \code{promax} from factor analysis. Hyperparameters for \code{varimax} and \code{promax} can be provided to \code{parameters}. Default is \code{TRUE}.
#' @param aggregate a logical value. If \code{TRUE}, bootstrap aggregation is applied to stabilize the method. Hyperparameters for bootstrap aggregation can be provided to \code{parameters}. Default is \code{TRUE}.
#' @param returnFull whether to return UHatList and VHatList. TODO!!!
#' @param parameters a list of hyperparameters.
#' \itemize{
#' \item \code{sprs.transform.m}: the power of the target used by \code{promax} if \code{sprs.transform} is \code{TRUE}. If \code{sprs.transform.m} is one, \code{varimax} is used. If \code{sprs.transform.m} is larger than one, \code{promax} is used with \code{sprs.transform.m}. Default is 1.
#' \item \code{U.prop.zero}: proportion of zeros in U's columns. For each column of U, we truncate the bottom \code{U.prop.zero} elements ranked by absolute values to zero. Default is 0.8.
#' \item \code{V.prop.zero}: proportion of zeros in V's columns. For each column of V, we truncate the bottom \code{V.prop.zero} elements ranked by absolute values to zero. Default is 0.8.
#' \item \code{aggregate.method}: method for aggregation if \code{aggregate} is \code{TRUE}. If \code{aggregate.method} is "bootstrap", SNPs are bootstrapped for aggregation. If \code{aggregate.method} is "subsample", SNPs are subsampled without replacement and the subsample size is determined by \code{subsample.p}. Default is "bootstrap".
#' \item \code{n.aggregate}: the number of aggregation samples. Default is 100.
#' \item \code{subsample.p}: subsample ratio if \code{aggregate.method} is "subsample". Default is 0.8.
#' \item \code{randomSeed}: random seed used in aggregation. Default is NULL.
#' \item \code{kU}: the number of genotypes in each genotype-phenotype cluster without aggregation. Default is the number of signal SNPs divided by 2\code{r}.
#' \item \code{kV}: the number of phenotypes in each genotype-phenotype cluster without aggregation. Default is the number of phenotypes divided by 2\code{r}.
#' \item \code{gg.weight}: weight of (genotype, genotype) pairs in computing co-appearance. Default is 1.
#' \item \code{gp.weight}: weight of (genotype, phenotype) pairs in computing co-appearance. Default is 1.
#' \item \code{pp.weight}: weight of (phenotype, phenotype) pairs in computing co-appearance. Default is 1.
#' \item \code{clustering.method}: method of clustering signal SNPs and phenotypes. If \code{aggregate} is \code{FALSE}, \code{clustering.method} is ignored and the i-th cluster consists of phenotypes and genotypes with non-zero elements in the i-th column of V and U, respectively.  If \code{aggregate} is \code{TRUE}, \code{clustering.method} currently takes the following options: if \code{clustering.method} equals "UMAP", we use \code{UMAP} to find two-dimensional embeddings of signal SNPs and use \code{kmeans} cluster the embeddings; if \code{clustering.method} equals "tSNE", we use \code{tSNE} to find two-dimensional embeddings of signal SNPs and phenotypes and use \code{kmeans} cluster the embeddings; if \code{clustering.method} equals "spectral", we use co-appearance as the kernel function and use spectral clustering to cluster signal SNPs and phenotypes. Hyperparameters for clustering can be further provided to \code{parameters}. Default is "UMAP".
#' \item \code{max.coappearance}: the upper bound of co-appearance weight. The distance induced by co-appearance is defined as max.coappearance - co-appearance. Default is the larger of 1/\code{r} and the maximal co-appearance.
#' \item \code{tSNE.perplexity}: TODO!!!
#' \item \code{UMAP.n_neighbors}: TODO!!!
#' \item \code{rsnps_results}: TODO!!!
#' \item \code{description}: TODO!!!
#' \item \code{pValue}: TODO!!!
#' \item \code{coappearance.truncate}: TODO!!!
#' }
#'
#' @return VHat, UHat
#'
#' @export
pathGPS = function(beta = "metabolites", UV = NULL, graph = NULL, p, p0, r, envr.subtract = TRUE, sprs.transform = TRUE, aggregate = TRUE, returnFull = FALSE,  parameters = list()){
  # preprocess
  if(is.character(beta)){
    if(beta == "metabolites"){
      beta = read.table(file = "./data/metabolitesBeta.txt", header = TRUE)
      rsnps_results = read.table("./data/metabolitesRsnps.txt", header = TRUE)
      description = read.table("./data/metabolitesDescription.txt", header = TRUE)
      suppressWarnings(RNGkind(sample.kind = "Rounding"))
    } else if(beta == "biobank"){
      beta = read.table(file = "./data/biobankBeta.txt", header = TRUE) # check.names = FALSE
      rsnps_results = read.table("./data/biobankRsnps.txt", header = TRUE)
      description = read.table("./data/biobankDescription.txt", header = TRUE)
      suppressWarnings(RNGkind(sample.kind = "Rounding"))
    }
  } else if(is.null(dim(beta))){stop("please input a matrix of summary statistics for multiple genotypes and multiple phenotypes")}
  pFull = dim(beta)[1]; q = dim(beta)[2]

  if(!is.null(UV)){
    UHatList = UV$UHat; VHatList = UV$VHat
    if(class(UHatList) != "list"){UHatList = list(UHatList)}
    if(class(VHatList) != "list"){VHatList = list(VHatList)}
    r = dim(UHatList[[1]])[2]
    print("Estimators of U, V are received. Next construct co-appearance graph...")
  } else {
    if(!is.null(colnames(beta))){phenotypeNames = colnames(beta)
    } else {phenotypeNames = paste("phenotype", seq(1,q), sep = "")}
    if(!is.null(rownames(beta))){SNPNames = rownames(beta)
    } else {SNPNames = paste("genotype", seq(1,pFull), sep = "")}

    if(sprs.transform){
      if(is.null(parameters$sprs.transform.m)){parameters$sprs.transform.m = 1}
      if(is.null(parameters$sprs.transform.lambda)){parameters$sprs.transform.lambda = 0}}
    if(is.null(parameters$V.prop.zero)){parameters$V.prop.zero = 0.8}
    if(is.null(parameters$U.prop.zero)){parameters$U.prop.zero = 0.8}
    if(!aggregate){parameters$n.aggregate = 1
    } else if(is.null(parameters$n.aggregate)){parameters$n.aggregate = 100}
    if(aggregate * (parameters$aggregate.method == "subsample")){
      if(is.null(parameters$pSample)){parameters$pSample = 2 * p}}
    if(!is.null(parameters$rsnps_results)){
      rsnps_results = parameters$rsnps_results; parameters$rsnps_results = NULL
    }
    if(!is.null(parameters$description)){
      description = parameters$description; parameters$description = NULL
    }
    if(!is.null(parameters$randomSeed)){set.seed(parameters$randomSeed)}

    # estimation
    VHatList = list(); UHatList = list()
    for(l in 1:parameters$n.aggregate){
      if (!aggregate){
        SNPIndex = seq(1, p) # no aggregation
      } else {
        if (parameters$aggregate.method == "subsample"){
          SNPIndex = sample(pSample, p, replace = FALSE)
        } else if (parameters$aggregate.method == "bootstrap"){
          SNPIndex = sample(p, p, replace = TRUE)
          uniqueSNPIndex = which(!duplicated(SNPIndex))
        }
      }
      beta1 = as.matrix(beta[SNPIndex, ], nrow = p) # beta is of dim at least 2x2
      beta0 = as.matrix(beta[seq(pFull + 1 - p0, pFull), ], nrow = p0)

      UVHat = pathGPS.oneshot(beta1 = beta1, beta0 = beta0, r = r, envr.subtract = envr.subtract, sprs.transform = sprs.transform, parameters = parameters)
      VHatList[[l]] = UVHat$VHat; UHatList[[l]] = UVHat$UHat
      if(parameters$aggregate.method == "bootstrap"){
        UHatList[[l]] = UHatList[[l]][uniqueSNPIndex,]
      }
    }
    print("Estimation of U, V done. Next construct co-appearance graph...")
  }

  # construct co-appearance graph
  if(is.null(parameters$kU)){parameters$kU = max(1, floor(p/2/r))}
  if(is.null(parameters$kV)){parameters$kV = max(1, floor(q/2/r))}
  if(is.null(parameters$gg.weight)){parameters$gg.weight = 1}
  if(is.null(parameters$gp.weight)){parameters$gp.weight = 1}
  if(is.null(parameters$pp.weight)){parameters$pp.weight = 1}
  g = EVToGraphFull(U = UHatList, V = VHatList, kU = parameters$kU, kV = parameters$kV, gp.weight = parameters$gp.weight, pp.weight = parameters$pp.weight, gg.weight = parameters$gg.weight)
  if(!is.null(parameters$coappearance.truncate)){
    weightMatrixZero = as.matrix(igraph::get.adjacency(g, type=c("both"), attr= "weight", names=TRUE))
    if(parameters$coappearance.truncate == TRUE){
      threshold = quantile(ecdf(weightMatrixZero), prob= 1 - 1/r)
    } else {threshold = parameters$coappearance.truncate}
    weightMatrixZero[weightMatrixZero <= threshold] = 0
    isolatedNode = colnames(weightMatrixZero)[which(apply(weightMatrixZero,2,sum) == 0)]
    if(length(isolatedNode) > 0) {g = igraph::delete.vertices(g, isolatedNode)}
  }
  print("Construction of co-appearance graph done. Next cluster genotypes and phenotypes...")

  # clustering
  weightMatrix = as.matrix(igraph::get.adjacency(g, type=c("both"), attr= "weight", names=TRUE)); diag(weightMatrix) = 0
  if(is.null(parameters$clustering.method)){parameters$clustering.method = "UMAP"}
  if(is.null(parameters$max.coappearance)){parameters$max.coappearance = max(1/r, max(weightMatrix))}
  if(!is.null(parameters$clustering.randomSeed)){set.seed(parameters$clustering.randomSeed)}
  if(parameters$clustering.method == "spectral"){
    kernMatrix = weightMatrix; diag(kernMatrix) = max(1/r, max(kernMatrix))
    clusteringMembership = myClustering(kernMatrix = kernMatrix, weightMatrix = weightMatrix, method = "spectral", parameters = parameters)$membership
    data = data.frame(x = rep(0, dim(kernMatrix)[1]),
                       y = rep(0, dim(kernMatrix)[1]),
                       group = as.factor(clusteringMembership),
                       name = colnames(kernMatrix))
  }
  if(parameters$clustering.method == "hierarchical"){
    distanceMatrix = parameters$max.coappearance - weightMatrix; diag(distanceMatrix) = 0
    clusteringMembership = myClustering(distMatrix = distanceMatrix, weightMatrix = weightMatrix, method = "hierarchical", parameters = parameters)$membership
    data = data.frame(x = rep(0, dim(kernMatrix)[1]),
                       y = rep(0, dim(kernMatrix)[1]),
                       group = as.factor(clusteringMembership),
                       name = colnames(kernMatrix))
  }
  if(parameters$clustering.method == "UMAP"){
    distanceMatrix = parameters$max.coappearance - weightMatrix; diag(distanceMatrix) = 0
    umap.parameters = umap::umap.defaults
    umap.parameters$n_components = 2; umap.parameters$input = "dist"
    if(!is.null(parameters$UMAP.n_neighbors)){
      umap.parameters$n_neighbors = parameters$UMAP.n_neighbors
    } else {umap.parameters$n_neighbors = floor(min(max(3, sqrt(dim(distanceMatrix)[1])), 15))}
    umapEmbedding = umap::umap(d = distanceMatrix, config = umap.parameters, method = "naive")
    clusteringMembership = myClustering(X = umapEmbedding$layout, weightMatrix = weightMatrix, method = "kmeans", parameters = parameters)$membership
    data = data.frame(x = umapEmbedding$layout[, 1],
                       y = umapEmbedding$layout[, 2],
                       group = as.factor(clusteringMembership),
                       name = colnames(distanceMatrix))
  }
  if(parameters$clustering.method == "tSNE"){
    distanceMatrix = parameters$max.coappearance - weightMatrix; diag(distanceMatrix) = 0
    if(is.null(parameters$tSNE.perplexity)){
      parameters$tSNE.perplexity = floor(min(max(5, sqrt(dim(distanceMatrix)[1])), (dim(distanceMatrix)[1]-1)/3))
    }
    distanceMatrix = parameters$max.coappearance - weightMatrix; diag(distanceMatrix) = 0
    tsne = Rtsne::Rtsne(X = distanceMatrix, is_distance = TRUE, perplexity = parameters$tSNE.perplexity)
    rownames(tsne$Y) = rownames(distanceMatrix)
    clusteringMembership = myClustering(X = tsne$Y, weightMatrix = weightMatrix, method = "kmeans", parameters = parameters)$membership
    data = data.frame(x = tsne$Y[, 1],
                       y = tsne$Y[, 2],
                       group = as.factor(clusteringMembership),
                       name = colnames(distanceMatrix))
  }
  data = data.table::data.table(data)
  data[, genotype := grepl("^rs", name)]

  if(!exists("rsnps_results")){print("no rsnps information");  rsnps_results = data.frame(data$name); colnames(description) = "rsid"; rsnps_results$gene = NA}
  if(!exists("description")){print("no phenotype description"); description = data.frame(data$name); colnames(description) = "name"; description$description = NA}

  data$position = seq(1, dim(weightMatrix)[1])
  data = merge(data, rsnps_results[, c("rsid", "gene")], by.x = "name", by.y = "rsid", all.x = TRUE)
  data[, plot.name := name]
  data[gene != "", plot.name := gene]
  colnames(weightMatrix)[data$position] = data$plot.name
  rownames(weightMatrix)[data$position] = data$plot.name
  clusteringSize = (table(clusteringMembership) + runif(table(clusteringMembership), 0, 0.01))[clusteringMembership]
  clusteringMembership = rep(seq(1,length(table(clusteringMembership))), sort(table(clusteringMembership)))
  names(clusteringMembership) = rownames(weightMatrix)[order(clusteringSize)]

  result = data[,c("name", "plot.name", "group", "genotype", "position", "x", "y")]
  colnames(result)[c(1,2,3)] = c("rsid", "name", "cluster")
  result$cluster[order(result$position)[order(clusteringSize)]] = rep(seq(1,length(table(clusteringMembership))), sort(table(clusteringMembership))) #???
  description$name[substr(description$name,1,1) %in% as.character(seq(0,9))] = paste("X", description$name[substr(description$name,1,1) %in% as.character(seq(0,9))], sep = "") #???

  result = merge(result, description, by.x = "name", by.y = "name", all.x = TRUE)
  result = merge(result, rsnps_results, by.x = "rsid", by.y = "rsid", all.x = TRUE)
  result = result[order(result$cluster),]
  result$description[result$genotype] = result$name[result$genotype]
  result$name[result$genotype] = result$rsid[result$genotype]
  result$rsid = NULL; result$gene = NULL; result$position = NULL
  print("Clustering of genotypes and phenotypes done.")

  result2 = list(); result2$cluster = result[,-c("x", "y")]; result2$graph = g; result2$embedding = data[,c("name", "x", "y")]; result2$parameters = parameters
  if(returnFull){result2$UHatList = UHatList; result2$VHatList = VHatList}
  return (result2)
}


#' pathGPS.oneshot
#'
#' oneshot estimation of U and V
#' TODO!!!
pathGPS.oneshot = function(beta1, beta0, r, envr.subtract = TRUE, sprs.transform = TRUE, parameters = list()){
  p = dim(beta1)[1]; q = dim(beta1)[2]; p0 = dim(beta0)[1]
  phenotypeNames = colnames(beta1); SNPNames = rownames(beta1)

  covBetaP = beta1 %*% t(beta1)
  covBetaQ1 = t(beta1) %*% beta1
  covBetaQ0 = t(beta0) %*% beta0

  # V
  if(envr.subtract){
    VHat= eigen(covBetaQ1 - covBetaQ0 * p / p0)$vectors[,seq(1,r)]
  } else {
    VHat = eigen(covBetaQ1)$vectors[,seq(1,r)]
  }

  # U
  UHat = beta1 %*% VHat %*% solve(t(VHat) %*% VHat)

  # linearly transform U, V to have sparser columns
  if(sprs.transform){
    if(parameters$sprs.transform.m == 1){
      UVTransform = biVarimax(U = UHat, V = VHat, lambda = parameters$sprs.transform.lambda, normalize = TRUE)
      VHat = UVTransform$loadingsV; UHat = UVTransform$loadingsU
    } else {
      UVTransform = biPromax(U = UHat, V = VHat, lambda = parameters$sprs.transform.lambda, normalize = TRUE, m = parameters$sprs.transform.m)
      VHat = UVTransform$loadingsV; UHat = UVTransform$loadingsU
    }
  }

  UHat = hardThresholding(UHat, zeroProp = parameters$U.prop.zero, normalize = FALSE)
  VHat = hardThresholding(VHat, zeroProp = parameters$V.prop.zero, normalize = FALSE)

  rownames(VHat) = phenotypeNames; rownames(UHat) = SNPNames
  result = list(); result$UHat = UHat; result$VHat = VHat
  result
}

