# helper of metabolites

#' colSpDistance
#'
#' compute the distance between the col space of X1 and that of X2
#' @param X covariate matrix to predict conditional densities at; each row represents an observation vector.
#'
#' @export
colSpDistance = function(X1,X2){
  if(is.null(dim(X1))){X1 = matrix(X1, ncol = 1); X2 = matrix(X2, ncol = 1)}
  if(is.null(dim(X2))){X1 = matrix(X1, ncol = 1); X2 = matrix(X2, ncol = 1)}
  if(rankMatrix(t(X1) %*% X1)[1] < dim(X1)[2]){dist = 0; return (dist)}
  if(rankMatrix(t(X2) %*% X2)[1] < dim(X2)[2]){dist = 0; return (dist)}
  P1 = X1 %*% solve(t(X1) %*% X1) %*% t(X1)
  P2 = X2 %*% solve(t(X2) %*% X2) %*% t(X2)
  dist = 1 - max(abs(eigen(P1-P2)$values))
  return (dist)
}

#' clusteringError
#'
#' compute the clustering error
#' @param V: phenotype - mediator matrix; each row represents a phenotype, each column represents a mediator. The number of mediators should be no larger than that of phenotypes.
#' @param VTrue: true phenotype - mediator matrix; each row represents a phenotype, each column represents a mediator. The number of mediators should be no larger than that of phenotypes.
#' @param U: genotype - mediator matrix estimator; each row represents a genotype, each column represents a mediator. The number of mediators should be no larger than that of genotypes. Default is NULL.
#' @param UTrue: true genotype - mediator matrix; each row represents a genotype, each column represents a mediator. The number of mediators should be no larger than that of genotypes. Default is NULL.
#'
#' @export
clusteringError = function(U = NULL, V, UTrue = NULL, VTrue, probUZero = 0, probVZero = 0){
  if(is.null(U)){
    # pre-process
    if(dim(V)[1] < dim(V)[2]){V = t(V)}
    if(dim(VTrue)[1] < dim(VTrue)[2]){VTrue = t(VTrue)}
    q = dim(V)[1]; r = max(dim(V)[2], dim(VTrue)[2]); rTrue = dim(VTrue)[2]
    if(is.null(probVZero)){probVZero = mean(VTrue == 0)}
    if(dim(V)[2] < dim(VTrue)[2]){
      V = cbind(V, matrix(0, nrow = q, ncol = dim(VTrue)[2] - dim(V)[2]))
    } else if(dim(V)[2] > dim(VTrue)[2]) {
      VTrue = cbind(VTrue, matrix(0, nrow = q, ncol = dim(V)[2] - dim(VTrue)[2]))
    }
    if(r > 6){stop("permutation is computationally prohibitive")}

    if(is.null(probVZero) == 0){
      V = hardThresholding(V, zeroProp = rep(probVZero, r), normalize = FALSE)
    }

    # all possible permutations
    permutation = rbind(seq(1,r), permute::allPerms(seq(1,r)))
    minError = q * r
    for(i in 1:dim(permutation)[1]){
      temp = clusteringErrorHelper(V = V[,permutation[i,]], VTrue = VTrue)
      if(minError > (temp$FP + temp$FN)){minError = temp$FP + temp$FN; maxPermute = permutation[i, ]}
    }

    # result
    result = clusteringErrorHelper(V = V[,maxPermute], VTrue = VTrue)
    result$error = minError/rTrue/q
    return(result)
  } else {
    # pre-process
    if(dim(U)[1] < dim(U)[2]){U = t(U)}
    if(dim(V)[1] < dim(V)[2]){V = t(V)}
    if(dim(UTrue)[1] < dim(UTrue)[2]){UTrue = t(UTrue)}
    if(dim(VTrue)[1] < dim(VTrue)[2]){VTrue = t(VTrue)}
    p = dim(U)[1]; r = max(dim(U)[2], dim(UTrue)[2]); q = dim(V)[1]; rTrue = dim(UTrue)[2]
    if(is.null(probUZero)){probUZero = mean(UTrue == 0)}
    if(is.null(probVZero)){probVZero = mean(VTrue == 0)}
    if(dim(V)[2] < dim(VTrue)[2]){
      U = cbind(U, matrix(0, nrow = p, ncol = dim(UTrue)[2] - dim(U)[2]))
      V = cbind(V, matrix(0, nrow = q, ncol = dim(VTrue)[2] - dim(V)[2]))
    } else if(dim(V)[2] > dim(VTrue)[2]) {
      UTrue = cbind(UTrue, matrix(0, nrow = p, ncol = dim(U)[2] - dim(UTrue)[2]))
      VTrue = cbind(VTrue, matrix(0, nrow = q, ncol = dim(V)[2] - dim(VTrue)[2]))
    }

    if(r > 6){stop("permutation is computationally prohibitive")}

    if((is.null(probUZero) + is.null(probVZero)) == 0){
      U = hardThresholding(U, zeroProp = rep(probUZero, r), normalize = FALSE)
      V = hardThresholding(V, zeroProp = rep(probVZero, r), normalize = FALSE)
    }

    # all possible permutations
    permutation = rbind(seq(1,r), permute::allPerms(seq(1,r)))
    minError = (p+q)*r
    for(i in 1:dim(permutation)[1]){
      temp = clusteringErrorHelper(U = U[,permutation[i,]], V = V[,permutation[i,]], UTrue = UTrue, VTrue = VTrue)
      if(minError > (temp$FP + temp$FN)){minError = temp$FP + temp$FN; maxPermute = permutation[i, ]}
    }

    # result
    result = clusteringErrorHelper(U = U[,maxPermute], V = V[,maxPermute], UTrue = UTrue, VTrue = VTrue)
    result$error = minError/rTrue/(p+q)
    return(result)
  }
}

clusteringErrorHelper = function(U = NULL, V, UTrue = NULL, VTrue){
  result = list()
  if(is.null(UTrue)){
    VNonZero = which(VTrue != 0)
    result$TP = sum(V[VNonZero] != 0)
    result$FP = sum(V[-VNonZero] != 0)
    result$TN = sum(V[-VNonZero] == 0)
    result$FN = sum(V[VNonZero] == 0)
  } else {
    UNonZero = which(UTrue != 0); VNonZero = which(VTrue != 0)
    result$TP = sum(U[UNonZero] != 0) + sum(V[VNonZero] != 0)
    result$FP = sum(U[-UNonZero] != 0) + sum(V[-VNonZero] != 0)
    result$TN = sum(U[-UNonZero] == 0) + sum(V[-VNonZero] == 0)
    result$FN = sum(U[UNonZero] == 0) + sum(V[VNonZero] == 0)
  }
  return(result)
}


# hard thresholding for sparse estimation
# input:
  # V: p by q vector, we truncate each column and preserve proportions of zeros in each column at zeroProp
  # zeroProp:proprotion of zeros for each column; if only a scalar is input, control the proportion of zeros in V
  # normalize: normalize the columns of the hard-thresholded V
hardThresholding = function(V, zeroProp, normalize = TRUE){
  if(is.null(V)){V = matrix(V, ncol = 1)}
  if(dim(V)[1] < dim(V)[2]){V = t(V)}
  p = dim(V)[1]; q = dim(V)[2]
  if(length(zeroProp) == 1){zeroProp = rep(zeroProp, q)}
  # if(length(zeroProp) == 1){
  #   index.max = apply(abs(V),2,which.max)
  #   index = order(abs(V), decreasing = TRUE)[1:floor((1-zeroProp) * p * q)]
  #   V[-index] = 0; V = matrix(V, ncol = q)
  #   if(normalize){
  #     VNorm = diag(t(V) %*% V)
  #     V[cbind(index.max[VNorm == 0], which(VNorm == 0))] = 1; VNorm[which(VNorm == 0)] = 1
  #     V = V %*% diag(1/sqrt(VNorm))
  #   }
  #   return(V)
  # }
  for(j in 1:q){
    nonZero.number = floor((1-zeroProp[j]) * p)
    if(nonZero.number <= 0){V[-which.max(abs(V[,j])), j] = 0}
    if(nonZero.number < p){
      nonZero.index = order(abs(V[,j]), decreasing = TRUE)[1:nonZero.number]
      V[-nonZero.index, j] = 0
    }
  }
  if(normalize){V = V %*% diag(1/sqrt(apply(V^2, 2, sum)))}
  return(V)
}

#' UVEstimator
#'
#' estimate the phenotype - mediator matrix V and the genotype - mediator matrix U
#' @param beta: non-null marginal association matrix; each row represents a non-null genotype, each column represents a phenotype
#' @param beta2: null marginal association matrix; each row represents a null genotype, each column represents a phenotype
#' @param r: number of latent mediators/pathways.
#' @param subtract: logical; if true, subtraction of null marginal association matrix from non-null marginal association matrix is conducted.
#' @param rotate: logical; if true, rotate the estimated U, V to have clear structures.
#' @param sparse: logical; if true, spca from package elasticnet is used to estimate sparse U, V.
#' @param hardThresholding: logical; if true, estimated U, V are hard-thresholded to preserve the desired number of non-zero entries.
#' @param parameters list of parameters.
#'
#' @export

UVEstimator = function(beta, beta2, r, subtract = TRUE, rotate = TRUE, hardThresholding = FALSE, sparse = FALSE, parameters = NULL){
  # pre-process
  p = dim(beta)[1]
  q = dim(beta)[2]
  p2 = dim(beta2)[1]

  covBetaP = beta %*% t(beta)
  covBetaQ = t(beta) %*% beta
  covBetaQ2 =  t(beta2) %*% beta2

  # V
  if(!sparse){
    if(!subtract){V = t(eigen(covBetaQ)$vectors[,seq(1,r)])
    } else {V = t(eigen(covBetaQ - covBetaQ2 * p/p2)$vectors[,seq(1,r)])}
    # U
    U = beta %*% t(V) %*% solve(V %*% t(V))
  } else {
    if(is.null(parameters$VVarnum)){stop("please input varnum for V")}
    if(is.null(parameters$UVarnum)){stop("please input varnum for U")}
    if(!subtract){
      V = t(elasticnet::spca(covBetaQ, K=r, para=rep(parameters$VVarnum, r), type=c("Gram"), sparse=c("varnum"), use.corr=FALSE, lambda=1e-6, max.iter=200, trace=FALSE, eps.conv=1e-3)$loadings)
    } else {
      V = t(elasticnet::spca(covBetaQ - covBetaQ2 * p/p2, K=r, para=rep(parameters$VVarnum, r), type=c("Gram"), sparse=c("varnum"), use.corr=FALSE, lambda=1e-6, max.iter=200, trace=FALSE, eps.conv=1e-3)$loadings)
    }
    U = beta %*% t(V) %*% solve(V %*% t(V))
    U = hardThresholding(U, zeroProp = rep(1-parameters$UVarnum/p, r), normalize = FALSE)
  }

  # rotate
  if(rotate){
    if(is.null(parameters)){
      parameters = list()
      parameters$rotateMethod = "varimax";
      parameters$rotateLambda = 0
      parameters$rotateM = 4
    }
    if(parameters$rotateMethod == "varimax"){
      tempRotate = biVarimax(U = U, V = V, lambda = parameters$rotateLambda, normalize = TRUE)
      V = t(tempRotate$loadingsV)
      U = tempRotate$loadingsU
    } else if (parameters$rotateMethod == "promax"){
      tempRotate = biPromax(U = U, V = V, lambda = parameters$rotateLambda, normalize = TRUE, m = parameters$rotateM)
      V = t(tempRotate$loadingsV)
      U = tempRotate$loadingsU
    }
  }

  if(hardThresholding){
    V = t(hardThresholding(V, zeroProp = 1-parameters$VVarnum/q, normalize = FALSE))
    U = hardThresholding(U, zeroProp = 1-parameters$UVarnum/p, normalize = FALSE)
  }

  # result
  result = list()
  result$U = U; result$V = V
  return(result)
}


#' biVarimax
#'
#' rotate a matrix pair (U, V) to have clear structures. Extended from varimax of a single matrix.
#' @param V: phenotype - mediator matrix; each row represents a phenotype, each column represents a mediator. The number of mediators should be no larger than that of phenotypes.
#' @param U: genotype - mediator matrix estimator; each row represents a genotype, each column represents a mediator. The number of mediators should be no larger than that of genotypes.
#' @param lambda hyper-parameter of balancing U and V: obj(U, V)  = lambda * obj_varimax(U) + obj_varimax(V). As lambda increases, more emphasis is put on U. Default is 1.
#' @param normalize logical; if true, the rows of U, V are re-scaled to unit length before rotation, and scaled back afterwards. Default is TRUE.
#' @param eps the tolerance for stopping: the relative change in the sum of singular values. Default is 1e-5.
#'
#' @export
biVarimax = function(U, V, lambda = 1, normalize = TRUE, eps = 1e-5){
  # preprocess
  if(is.null(dim(U)) || is.null(dim(V))){stop("please input a matrix pair")}
  if(dim(U)[1] < dim(U)[2]){U = t(U)}
  if(dim(V)[1] < dim(V)[2]){V = t(V)}
  r = dim(U)[2]; p = dim(U)[1]; q = dim(V)[1]
  if(r != dim(V)[2]){stop("please input a matrix pair of the right dimension")}

  # normalize
  if(normalize){
    scU = sqrt(drop(apply(U, 1L, function(x) sum(x^2)))); scU[which(scU == 0)] = 1; U = U/scU
    scV = sqrt(drop(apply(V, 1L, function(x) sum(x^2)))); scV[which(scV == 0)] = 1; V = V/scV
  }

  # BSV algorithm with alpha = 0
  R = diag(r); d = 0
  for (i in 1L:1000L){
    UR = U %*% R
    VR = V %*% R
    BU = t(U) %*% (UR^3 - UR %*% diag(drop(rep(1, p) %*% UR^2))/p)/p
    BV = t(V) %*% (VR^3 - VR %*% diag(drop(rep(1, q) %*% VR^2))/q)/q
    B = lambda * BU + BV
    sB = La.svd(B)
    R = sB$u %*% sB$vt
    dNew = sum(sB$d)
    if (dNew < d * (1 + eps)){break} else {d = dNew}
  }
  UR = U %*% R
  VR = V %*% R

  # post-process
  if(normalize){UR = UR * scU; VR = VR * scV}
  dimnames(UR) <- dimnames(U); dimnames(VR) <- dimnames(V)
  result = list()
  result$loadingsU = UR; result$loadingsV = VR; result$rotmat = R
  return(result)
}


#' biPromax
#'
#' rotate a matrix pair (U, V) to have clear structures. Extended from promax of a single matrix.
#' @param V: phenotype - mediator matrix; each row represents a phenotype, each column represents a mediator. The number of mediators should be no larger than that of phenotypes.
#' @param U: genotype - mediator matrix estimator; each row represents a genotype, each column represents a mediator. The number of mediators should be no larger than that of genotypes.
#' @param lambda hyper-parameter of balancing U and V: obj(U, V)  = lambda * obj_varimax(U) + obj_varimax(V). As lambda increases, more emphasis is put on U. Default is 1.
#' @param normalizeVarimax logical; if true, the rows of U, V are re-scaled to unit length before rotation, and scaled back afterwards. Default is TRUE.
#' @param eps the tolerance for stopping: the relative change in the sum of singular values. Default is 1e-5.
#' @param m the power used the target for promax. Values of 2 to 4 are recommended. Default is 4.
#'
#' @export
biPromax = function(U, V, lambda = 1, normalizeVarimax = TRUE, eps = 1e-05, m = 4){
  # preprocess
  if(is.null(dim(U)) || is.null(dim(V))){stop("please input a matrix pair")}
  if(dim(U)[1] < dim(U)[2]){U = t(U)}
  if(dim(V)[1] < dim(V)[2]){V = t(V)}
  r = dim(U)[2]; p = dim(U)[1]; q = dim(V)[1]
  if(r != dim(V)[2]){stop("please input a matrix pair with matching smaller dimensions")}

  # bi-Varimax
  temp = biVarimax(U = U, V = V, lambda = lambda, normalize = normalizeVarimax, eps = eps)

  # bi-Promax
  UR = temp$loadingsU
  VR = temp$loadingsV
  QR = UR * abs(UR)^(m-1)
  PR = VR * abs(VR)^(m-1)

  # initialize
  lambdaLC = 1 # hyper-parameter for ||ML^\top - I||_F^2
  M = diag(rep(1,r)); L = diag(rep(1,r))
  # alternating
  for (i in 1:10){
    for(j in 1:100){
      L = solve(t(VR) %*% VR + lambdaLC * M %*% t(M)) %*% (t(VR) %*% PR + lambdaLC * M)
      M = solve(lambda * t(UR) %*% UR + lambdaLC * L %*% t(L)) %*% (lambda * t(UR) %*% QR + lambdaLC * L)
    }
    lambdaLC = 2 * lambdaLC
  }
  # normalize
  d = diag(solve(t(L) %*% L))
  L = temp$rotmat %*% L %*% diag(sqrt(d))
  M = temp$rotmat %*% M %*% diag(1/sqrt(d))
  VP  = V %*% L
  UP  = U %*% M

  dimnames(UP) <- dimnames(U); dimnames(VP) <- dimnames(V)
  result = list()
  result$loadingsU = UP; result$loadingsV = VP; result$rotmatU = M; result$rotmatV = L
  return(result)
}

#' # transfrorm EVs to a weighted undirected graph
#' #' @export
#' EVToGraph = function(V, varProp = NULL, k = NULL, eigenvalues = NULL, returnFull = FALSE){
#'   edgeList = list()
#'   if(class(V) != "list"){
#'     r = min(dim(V))
#'     edgeList = EVToGraphHelper(V = V, edgeList = edgeList,
#'                                varProp = varProp, k = k, eigenvalues = eigenvalues)
#'   } else {
#'     r = min(dim(V[[1]]))
#'     for(i in 1:length(V)){
#'       edgeList = EVToGraphHelper(V = V[[i]], edgeList = edgeList,
#'                                  varProp = varProp, k = k, eigenvalues = eigenvalues)
#'     }
#'   }
#'   # construct graphs
#'   if(length(edgeList) == 0){print("no edge"); return(NULL)}
#'   weight = table(unlist(edgeList))/r; if(class(V) == "list"){weight = weight/length(V)}
#'   edge = unlist(strsplit(names(weight),split = ";")); g = igraph::graph(edge) # isolates not needed
#'   names(weight) = NULL; igraph::E(g)$weight = weight
#'   g = igraph::as.undirected(g, mode= "collapse",edge.attr.comb=list(weight="sum", "ignore"))
#'   if(!returnFull){return(g)
#'     } else {result = list(); result$g = g; result$edge = edge; result$weight = weight; return(result)}
#' }

# EVToGraphHelper = function(V, edgeList, varProp = NULL, k = NULL, eigenvalues = NULL){
#   # preprocess
#   if(is.null(dim(V))){V = matrix(V, ncol = 1)}
#   if(dim(V)[1] <= dim(V)[2]){V = t(V)}
#   q = dim(V)[1]; r = dim(V)[2]
#   if(is.null(rownames(V))){stop("please input names of V")} else {nodeNames = rownames(V)}
#   edgeCounter = length(edgeList)
#   if(is.null(eigenvalues)){
#     # determine dominating phenotypes
#     if(!is.null(k)){
#       for(i in 1:r){
#         index = order((V[,i])^2, decreasing = TRUE)[1:k]
#         for(j in 2:k){
#           for(l in 1:(j-1)){
#             edgeCounter = edgeCounter + 1
#             edgeList[[edgeCounter]] = paste(nodeNames[index[j]], nodeNames[index[l]], sep = ";")
#           }
#         }
#       }
#     } else if(!is.null(varProp)){
#       for(i in 1:r){
#         maxIndex = suppressWarnings(max(which(cumsum(sort((V[,i])^2, decreasing = TRUE)) <= varProp)))
#         if(maxIndex <= 1){next}
#         index = order((V[,i])^2, decreasing = TRUE)[1:maxIndex]
#         for(j in 2:maxIndex){
#           for(l in 1:(j-1)){
#             edgeCounter = edgeCounter + 1
#             edgeList[[edgeCounter]] = paste(nodeNames[index[j]], nodeNames[index[l]], sep = ";")
#           }
#         }
#       }
#     }
#   }
#   return (edgeList)
# }

#' #' @export
#' weightToGraph = function(weightMatrix, mode = "undirected"){
#'   diag(weightMatrix) = 0
#'   g = igraph::graph_from_adjacency_matrix(weightMatrix, mode = mode, weighted = TRUE)
#'   g = igraph::as.undirected(g, mode= "collapse",edge.attr.comb=list(weight="sum", "ignore"))
#'   return(g)
#' }

# transfrorm EVs to a full weighted undirected graph with three types of edges
#' @export
EVToGraphFull = function(U, V, kU = NULL,kV = NULL, gp.weight = 1, pp.weight = 1, gg.weight = 1, returnFull = FALSE, weightMatrixBootstrap = NULL){
  edgeList = list(); edgeList$edgeListUU = list(); edgeList$edgeListVV = list(); edgeList$edgeListUV = list()
  if(class(V) != "list"){V = list(V)}
  if(class(U) != "list"){U = list(U)}
  # if(class(V) != "list"){
  #   r = min(dim(V))
  #   edgeList = EVToGraphFullHelper(U = U, V = V, edgeList = edgeList,
  #                              varPropU = varPropU, kU = kU, varPropV = varPropV,
  #                              kV = kV, eigenvalues = eigenvalues)
  # } else {
  #   r = min(dim(V[[1]]))
  #   for(i in 1:length(V)){
  #     edgeList = EVToGraphFullHelper(U = U[[i]], V = V[[i]], edgeList = edgeList,
  #                                   varPropU = varPropU, kU = kU, varPropV =
  #                                   varPropV, kV = kV, eigenvalues = eigenvalues)
  #   }
  # }
  r = min(dim(V[[1]]))
  for(i in 1:length(V)){
    edgeList = EVToGraphFullHelper(U = U[[i]], V = V[[i]], edgeList = edgeList, kU = kU, kV = kV)
  }

  if((length(edgeList$edgeListUU) + length(edgeList$edgeListVV) + length(edgeList$edgeListUV)) == 0){print("no edge detected"); return(NULL)}
  weight = c(table(unlist(edgeList$edgeListUU)) * gg.weight, table(unlist(edgeList$edgeListVV)) * pp.weight, table(unlist(edgeList$edgeListUV)) * gp.weight)
  indexZeroWeight = which(weight == 0)
  if(length(indexZeroWeight) != 0){weight = weight[-indexZeroWeight]}
  edge = unlist(strsplit(names(weight),split = ";")); g = igraph::graph(edge)
  weight = weight/r/length(V); names(weight) = NULL; igraph::E(g)$weight = weight
  g = igraph::as.undirected(g, mode= "collapse", edge.attr.comb=list(weight="sum", "ignore"))

  # remove non-significant edges
  # if(!is.null(weightMatrixBootstrap)){
  #   weightMatrix = igraph::get.adjacency(g, type=c("both"), attr="weight", names=TRUE); weightMatrix = as.matrix(weightMatrix)
  #   weightMatrixZero = weightMatrix * weightMatrixBootstrap[colnames(weightMatrix), colnames(weightMatrix)]
  #   isolatedNode = colnames(weightMatrix)[which(apply(weightMatrixZero,2,sum) == 0)]
  #   if(length(isolatedNode) > 0) {g = igraph::delete.vertices(g, isolatedNode)}
  # }

  if(!returnFull){return(g)
  } else {
    result = list(); result$g = g; result$edge = edge; result$weight = weight
    return(result)
  }
}

EVToGraphFullHelper = function(U, V, edgeList, kU = NULL, kV = NULL){
  # preprocess
  if(is.null(dim(V))){V = matrix(V, ncol = 1)}
  if(is.null(dim(V))){U = matrix(U, ncol = 1)}
  if(dim(V)[1] <= dim(V)[2]){V = t(V)}
  if(dim(U)[1] <= dim(U)[2]){U = t(U)}
  q = dim(V)[1]; r = dim(V)[2]; p = dim(U)[1]
  if(is.null(rownames(V))){stop("please input rownames of V")} else {nodeNamesV = rownames(V)}
  if(is.null(rownames(U))){stop("please input rownames of U")} else {nodeNamesU = rownames(U)}
  edgeCounterVV = length(edgeList$edgeListVV)
  edgeCounterUV = length(edgeList$edgeListUV)
  edgeCounterUU = length(edgeList$edgeListUU)

  # determine dominating phenotypes
  for(i in 1:r){
    indexU = order((U[,i])^2, decreasing = TRUE)[1:min(kU, sum(U[,i] != 0))]
    indexV = order((V[,i])^2, decreasing = TRUE)[1:min(kV, sum(V[,i] != 0))]
    # VV
    if(length(indexV) > 1){
      for(j in 2:length(indexV)){
        for(l in 1:(j-1)){
          edgeCounterVV = edgeCounterVV + 1
          edgeList$edgeListVV[[edgeCounterVV]] = paste(nodeNamesV[indexV[j]], nodeNamesV[indexV[l]], sep = ";")
        }
      }
    }
    # UU
    if(length(indexU) > 1){
      for(j in 2:length(indexU)){
        for(l in 1:(j-1)){
          edgeCounterUU = edgeCounterUU + 1
          edgeList$edgeListUU[[edgeCounterUU]] = paste(nodeNamesU[indexU[j]], nodeNamesU[indexU[l]], sep = ";")
        }
      }
    }
    # UV
    for(j in 1:length(indexU)){
      for(l in 1:length(indexV)){
        edgeCounterUV = edgeCounterUV + 1
        edgeList$edgeListUV[[edgeCounterUV]] = paste(nodeNamesU[indexU[j]], nodeNamesV[indexV[l]], sep = ";")
      }
    }
  }
  # # UU
  # for(i in 1:r){
  #   indexU = order((U[,i])^2, decreasing = TRUE)[1:min(kU, sum(U[,i]^2 != 0))]
  #   if(length(indexU) > 1){
  #     for(j in 2:length(indexU)){
  #       for(l in 1:(j-1)){
  #         edgeCounterUU = edgeCounterUU + 1
  #         edgeList$edgeListUU[[edgeCounterUU]] = paste(nodeNamesU[indexU[j]], nodeNamesU[indexU[l]], sep = ";")
  #       }
  #     }
  #   }
  # }
  # # UV
  # for(i in 1:r){
  #   indexU = order((U[,i])^2, decreasing = TRUE)[1:min(kU, sum(U[,i]^2 != 0))]
  #   indexV = order((V[,i])^2, decreasing = TRUE)[1:min(kV, sum(V[,i]^2 != 0))]
  #   for(j in 1:length(indexU)){
  #     for(l in 1:length(indexV)){
  #       edgeCounterUV = edgeCounterUV + 1
  #       edgeList$edgeListUV[[edgeCounterUV]] = paste(nodeNamesU[indexU[j]], nodeNamesV[indexV[l]], sep = ";")
  #     }
  #   }
  # }
  return (edgeList)
}


# bipartite graph
# EVToGraphBipartite = function(V, varProp = NULL, k = NULL, eigenvalues = NULL, returnFull = FALSE){
#   if(class(V) != "list"){
#     edgeList = EVToGraphHelper(V = V, edgeList = list(),
#                                varProp = varProp, k = k, eigenvalues = eigenvalues)
#   } else {
#     edgeList = list()
#     for(i in 1:length(V)){
#       edgeList = EVToGraphHelper(V = V[[i]], edgeList = edgeList,
#                                  varProp = varProp, k = k, eigenvalues = eigenvalues)
#     }
#   }
#   # construct graphs
#   if(length(edgeList) == 0){print("no edge"); return(NULL)}
#   weight = table(unlist(edgeList))/r; if(class(V) == "list"){weight = weight/length(V)}
#   edge = unlist(strsplit(names(weight),split = ";")); g = igraph::graph(edge) # isolates not needed
#   names(weight) = NULL; igraph::E(g)$weight = weight
#   g = igraph::as.undirected(g, mode= "collapse",edge.attr.comb=list(weight="sum", "ignore"))
#   if(!returnFull){return(g)
#   } else {result = list(); result$g = g; result$edge = edge; result$weight = weight; return(result)}
# }

# EVToGraphBipartiteHelper = function(V, edgeList, varProp = NULL, k = NULL, eigenvalues = NULL){
#   # preprocess
#   if(is.null(dim(V))){V = matrix(V, ncol = 1)}
#   if(dim(V)[1] <= dim(V)[2]){V = t(V)}
#   q = dim(V)[1]; r = dim(V)[2]
#   if(is.null(rownames(V))){stop("please input names of V")} else {nodeNames = rownames(V)}
#   edgeCounter = length(edgeList)
#   if(is.null(eigenvalues)){
#     # determine dominating phenotypes
#     if(!is.null(k)){
#       for(i in 1:r){
#         index = order((V[,i])^2, decreasing = TRUE)[1:k]
#         for(j in 2:k){
#           for(l in 1:(j-1)){
#             edgeCounter = edgeCounter + 1
#             edgeList[[edgeCounter]] = paste(nodeNames[index[j]], nodeNames[index[l]], sep = ";")
#           }
#         }
#       }
#     } else if(!is.null(varProp)){
#       for(i in 1:r){
#         maxIndex = suppressWarnings(max(which(cumsum(sort((V[,i])^2, decreasing = TRUE)) <= varProp)))
#         if(maxIndex <= 1){return (edgeList)}
#         index = order((V[,i])^2, decreasing = TRUE)[1:maxIndex]
#         for(j in 2:maxIndex){
#           for(l in 1:(j-1)){
#             edgeCounter = edgeCounter + 1
#             edgeList[[edgeCounter]] = paste(nodeNames[index[j]], nodeNames[index[l]], sep = ";")
#           }
#         }
#       }
#     }
#   }
#   return (edgeList)
# }


# clustering methods with hyper-parameter tuning
# input:
#   X: named matrix of features: each row represents an element, each column represents a feature
#   distMatrix: named pairwise distance matrix
#   kernMatrix: named kernel matrix
#   method: clustering method: kmeans, cmeans, spectral, hierarchical
#   parameters: list of candidate tuning parameters of the the clustering method: gg.weight, gp.weight, pp.weight, similarity.distance, number.center
# output:
#   weightMatrix: named weightMatrix to compare to in tuning
#   membership: matrix of membership: each row represents a cluster, each column represents a named element
myClustering = function(X = NULL, distMatrix = NULL, kernMatrix = NULL, weightMatrix, method, parameters = list()){
  # preprocess
  if(is.null(parameters$similarity.distance)){parameters$similarity.distance = "Frobenius"}

  SNPIndex = which(grepl("^rs", rownames(weightMatrix)))
  phenotypeIndex = which(!grepl("^rs", rownames(weightMatrix)))

  # clustering
  if(method == "kmeans"){
    if(!is.null(parameters$number.center)){
      numberCenterSeq = floor(pmax(2, pmin(parameters$number.center,  dim(X)[1]/2-1)))} else {
      numberCenterSeq = floor(pmax(2, pmin(seq(2, sqrt(dim(X)[1]), length.out = 5), dim(X)[1]/2-1)))
    }
    similarityDistance = rep(0, length(numberCenterSeq))
    for(i in 1:length(numberCenterSeq)){
      if(!is.null(parameters$randomSeed2)){set.seed(parameters$randomSeed2)}
      kmeansMembership = kmeans(X, centers = numberCenterSeq[i], nstart = 20)$cluster
      weightMatrixTemp = matrix(rep(kmeansMembership, length(kmeansMembership)), length(kmeansMembership))
      weightMatrixTemp = (weightMatrixTemp == t(weightMatrixTemp))
      diag(weightMatrixTemp) = 0
      if(length(SNPIndex) > 0){
        weightMatrixTemp[SNPIndex, SNPIndex] = weightMatrixTemp[SNPIndex, SNPIndex] * parameters$gg.weight
        if(length(phenotypeIndex) > 0){
          weightMatrixTemp[SNPIndex, phenotypeIndex] = weightMatrixTemp[SNPIndex, phenotypeIndex] * parameters$gp.weight
          weightMatrixTemp[phenotypeIndex, SNPIndex] = weightMatrixTemp[phenotypeIndex, SNPIndex] * parameters$gp.weight
        }
      }
      if(length(phenotypeIndex) > 0){
        weightMatrixTemp[phenotypeIndex, phenotypeIndex] = weightMatrixTemp[phenotypeIndex, phenotypeIndex] * parameters$pp.weight
      }
      if(parameters$similarity.distance == "Frobenius"){
        lmTemp = lm(c(weightMatrix) ~ c(weightMatrixTemp) - 1)
        similarityDistance[i] = 1-summary(lmTemp)[[8]]
      } else if(parameters$similarity.distance == "hamming"){
        similarityDistance[i] = mean((weightMatrix!=0) != (weightMatrixTemp!=0))
      }
    }
    index = which.min(similarityDistance)
    if(!is.null(parameters$randomSeed2)){set.seed(parameters$randomSeed2)}
    membership = kmeans(X, centers = numberCenterSeq[index], nstart = 20)$cluster
    V = matrix(0, nrow=numberCenterSeq[index], ncol = length(membership))
    colnames(V) = names(membership)
    V[cbind(membership, seq(1, length(membership)))] = 1
  } else if(method == "spectral"){
    temp = kernlab::kernelMatrix(kernlab::rbfdot(), matrix(rnorm((dim(kernMatrix)[1])^2),dim(kernMatrix)[1], dim(kernMatrix)[1]))
    temp[temp != Inf] = kernMatrix
    if(!is.null(parameters$number.center)){
      numberCenterSeq = floor(pmax(3, pmin(parameters$number.center, dim(kernMatrix)[1]/2-1)))
    } else {
      numberCenterSeq = floor(pmax(3, pmin(seq(3, sqrt(dim(kernMatrix)[1]),length.out = 5), dim(kernMatrix)[1]/2-1)))
    }
    similarityDistance = rep(0, length(numberCenterSeq))
    for(i in 1:length(numberCenterSeq)){
      if(!is.null(parameters$randomSeed2)){set.seed(parameters$randomSeed2)}
      spcMembership = kernlab::specc(x = temp, centers = numberCenterSeq[i])
      weightMatrixTemp = matrix(rep(spcMembership, length(spcMembership)), length(spcMembership))
      weightMatrixTemp = (weightMatrixTemp == t(weightMatrixTemp))
      diag(weightMatrixTemp) = 0
      if(length(SNPIndex) > 0){
        weightMatrixTemp[SNPIndex, SNPIndex] = weightMatrixTemp[SNPIndex, SNPIndex] * parameters$gg.weight
        if(length(phenotypeIndex) > 0){
          weightMatrixTemp[SNPIndex, phenotypeIndex] = weightMatrixTemp[SNPIndex, phenotypeIndex] * parameters$gp.weight
          weightMatrixTemp[phenotypeIndex, SNPIndex] = weightMatrixTemp[phenotypeIndex, SNPIndex] * parameters$gp.weight
        }
      }
      if(length(phenotypeIndex) > 0){
        weightMatrixTemp[phenotypeIndex, phenotypeIndex] = weightMatrixTemp[phenotypeIndex, phenotypeIndex] * parameters$pp.weight
      }
      if(parameters$similarity.distance == "Frobenius"){
        lmTemp = lm(c(weightMatrix) ~ c(weightMatrixTemp) - 1)
        similarityDistance[i] = 1-summary(lmTemp)[[8]]
      } else if(parameters$similarity.distance == "hamming"){
        similarityDistance[i] = mean((weightMatrix!=0) != (weightMatrixTemp!=0))
      }
    }
    index = which.min(similarityDistance)
    if(!is.null(parameters$randomSeed2)){set.seed(parameters$randomSeed2)}
    membership = kernlab::specc(x = temp, centers = numberCenterSeq[index])
    names(membership) = rownames(kernMatrix)
    V = matrix(0, nrow=numberCenterSeq[index], ncol = length(membership))
    colnames(V) = names(membership)
    V[cbind(membership, seq(1, length(membership)))] = 1
  } else if(method == "hierarchical"){
    hier = hclust(as.dist(distMatrix), method = "complete")
    if(is.null(parameters$number.center)){
      numberCenterSeq = floor(pmax(2, pmin(parameters$number.center,  dim(distMatrix)[1]/2-1)))
    } else{
      numberCenterSeq = floor(pmax(2, pmin(seq(2, sqrt(dim(distMatrix)[1]), length.out = 5),  dim(distMatrix)[1]/2-1)))}
    similarityDistance = rep(0, length(numberCenterSeq))
    for(i in 1:length(numberCenterSeq)){
      hierMembership = cutree(hier, k = numberCenterSeq[i])
      weightMatrixTemp = matrix(rep(hierMembership, length(hierMembership)), length(hierMembership))
      weightMatrixTemp = (weightMatrixTemp == t(weightMatrixTemp))
      diag(weightMatrixTemp) = 0
      if(length(SNPIndex) > 0){
        weightMatrixTemp[SNPIndex, SNPIndex] = weightMatrixTemp[SNPIndex, SNPIndex] * parameters$gg.weight
        if(length(phenotypeIndex) > 0){
          weightMatrixTemp[SNPIndex, phenotypeIndex] = weightMatrixTemp[SNPIndex, phenotypeIndex] * parameters$gp.weight
          weightMatrixTemp[phenotypeIndex, SNPIndex] = weightMatrixTemp[phenotypeIndex, SNPIndex] * parameters$gp.weight
        }
      }
      if(length(phenotypeIndex) > 0){
        weightMatrixTemp[phenotypeIndex, phenotypeIndex] = weightMatrixTemp[phenotypeIndex, phenotypeIndex] * parameters$pp.weight
      }
      if(parameters$similarity.distance == "Frobenius"){
        lmTemp = lm(c(weightMatrix) ~ c(weightMatrixTemp) - 1)
        similarityDistance[i] = 1-summary(lmTemp)[[8]]
      } else if(parameters$similarity.distance == "hamming"){
        similarityDistance[i] = mean((weightMatrix!=0) != (weightMatrixTemp!=0))
      }
    }
    index = which.min(similarityDistance)
    membership = cutree(hier, k = numberCenterSeq[index])
    names(membership) = rownames(distMatrix)
    V = matrix(0, nrow=numberCenterSeq[index], ncol = length(membership))
    colnames(V) = names(membership)
    V[cbind(membership, seq(1, length(membership)))] = 1
  } else if(method == "cmeans"){
    if(!is.null(parameters$number.center)){
      numberCenterSeq = floor(pmax(2, pmin(parameters$number.center,  dim(X)[1]/2-1)))} else {
      numberCenterSeq = floor(pmax(2, pmin(seq(2, sqrt(dim(X)[1]), length.out = 5),  dim(X)[1]/2-1)))}
    if(!is.null(parameters$memb.exp)){
      memb.expSeq = floor(pmax(2, pmin(parameters$memb.exp,  dim(X)[1]/2-1)))} else {memb.expSeq = c(1.5,2)}
    if(!is.null(parameters$threshold)){
      thresholdSeq = floor(pmax(2, pmin(parameters$thresholdSeq,  dim(X)[1]/2-1)))} else {thresholdSeq = 1/2^(seq(2,6))}
    similarityDistance = array(0, dim = c(length(numberCenterSeq), length(memb.expSeq), length(thresholdSeq)))
    for(i in 1:length(numberCenterSeq)){
      for(j in 1:length(memb.expSeq)){
        probMatrix = cluster::fanny(x = X, k = numberCenterSeq[i], diss = FALSE, memb.exp= memb.expSeq[j])
        for(l in 1:length(thresholdSeq)){
          probMatrixTemp = (probMatrix$membership >= thresholdSeq[l])
          weightMatrixTemp = probMatrixTemp %*% t(probMatrixTemp)
          diag(weightMatrixTemp) = 0
          if(length(SNPIndex) > 0){
            weightMatrixTemp[SNPIndex, SNPIndex] = weightMatrixTemp[SNPIndex, SNPIndex] * parameters$gg.weight
            if(length(phenotypeIndex) > 0){
              weightMatrixTemp[SNPIndex, phenotypeIndex] = weightMatrixTemp[SNPIndex, phenotypeIndex] * parameters$gp.weight
              weightMatrixTemp[phenotypeIndex, SNPIndex] = weightMatrixTemp[phenotypeIndex, SNPIndex] * parameters$gp.weight
            }
          }
          if(length(phenotypeIndex) > 0){
            weightMatrixTemp[phenotypeIndex, phenotypeIndex] = weightMatrixTemp[phenotypeIndex, phenotypeIndex] * parameters$pp.weight
          }
          if(parameters$similarity.distance == "Frobenius"){
            lmTemp = lm(c(weightMatrix) ~ c(weightMatrixTemp) - 1)
            similarityDistance[i,j,l] = 1-summary(lmTemp)[[8]]
          } else if(parameters$similarity.distance == "hamming"){
            similarityDistance[i,j,l] = mean((weightMatrix!=0) != (weightMatrixTemp!=0))
          }
        }
      }
    }
    index = which(similarityDistance == min(similarityDistance), arr.ind = TRUE)[1,]
    membership = cluster::fanny(x = X, k = numberCenterSeq[index[1]], diss = FALSE, memb.exp= memb.expSeq[index[2]])$membership
    V = t(membership >= thresholdSeq[index[3]])
    colnames(V) = rownames(X)
  }
  result = list(); result$V = V; result$similarity.distance.seq = similarityDistance; result$membership = membership
  return (result)
}

# refFile: .tsv., header =TRUE. 1: ID; 2: description.
readDescription = function(refFile, membership = NULL, names = NULL){
  test1 = read.table(refFile, sep = "\t", header = TRUE)[,1:2]
  colnames(test1) = c("name", "description")
  if(is.null(names)){names = names(membership)}
  names = createNewNames(names)
  test2 = as.data.frame(as.numeric(membership)); colnames(test2) = "membership";  test2$name = names; rownames(test2) = NULL
  test = merge(test1, test2, all.x = TRUE, by.x = "name", by.y = "name")
  return(test)
}

createNewNames = function(names){
  for (i in 1:length(names)){
    if(substr(names[i],1,1) == "X"){
      if(substr(names[i],1,2) != "XI"){
        names[i] = substr(names[i],2,1000)
      }
    }
  }
  return(names)
}

