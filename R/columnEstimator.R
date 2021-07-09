# helper of metabolites

#' columnEstimator
#'
#' estimate the row space of V (#mediators x #phenotypes) and the column space of U (#genotypes x #mediators)
#' @param X covariate matrix to predict conditional densities at; each row represents an observation vector.
#' @param betas
#' @param p
#' @param pNull
#' @param pSample
#' @param r
#' @param V.subtract
#' @param rotate
#' @param U.project
#' @param UV.truncate
#' @param aggregate.method: "subsample", "bootstrap"
#' @param parameters: list of parameters: rotate.m, propUZero, probVZero, randomSeed, similarity.truncate, n.subsample
#'
#' @return VHat, UHat
#'
#' @export
columnEstimator = function(betas = "metabolites", p, pNull, pSample = NULL, r, V.subtract = TRUE, rotate = TRUE, rotate.m = 1, U.project = TRUE, UV.truncate = FALSE, aggregate.method = "bootstrap", parameters = list()){
  # preprocess
  if(is.character(betas)){
    if(betas == "metabolites"){
      betas = read.table(file = "./data/metabolitesBeta.txt", header = TRUE)
      suppressWarnings(RNGkind(sample.kind = "Rounding"))
    } else if(betas == "biobankSub"){
      betas = read.table(file = "./data/biobankSubBeta.txt", header = TRUE) # check.names = FALSE
      suppressWarnings(RNGkind(sample.kind = "Rounding"))
    } else if(betas == "biobank"){
      betas = read.table(file = "./data/biobankBeta.txt", header = TRUE) # check.names = FALSE
      suppressWarnings(RNGkind(sample.kind = "Rounding"))
    }
  } else if(is.null(dim(betas))){stop("please input a matrix of summary statistics for multiple genotypes and multiple phenotypes")}
  pTotal = dim(betas)[1]; q = dim(betas)[2]
  if(is.null(pSample)){pSample = 2 * p}

  if(!is.null(colnames(betas))){phenotypeNames = colnames(betas)
  } else {phenotypeNames = paste("phenotype", seq(1,q), sep = "")}
  if(!is.null(rownames(betas))){SNPNames = rownames(betas)
  } else {SNPNames = paste("genotype", seq(1,pTotal), sep = "")}

  if(is.null(parameters$n.subsample)){parameters$n.subsample = 100}
  if(is.null(parameters$similarity.truncate)){parameters$similarity.truncate = FALSE} else if(parameters$similarity.truncate){
    if(aggregate.method == "subsample"){stop("similarity truncation is only available for bootstrapping")}
    resultDist = array(0, dim = c(parameters$n.subsample, q + p, q + p))
  }

  if(is.null(parameters$randomSeed)){parameters$randomSeed = 318}
  set.seed(parameters$randomSeed)

  # estimation
  VHatList = list(); UHatList = list()
  for(l in 1:parameters$n.subsample){
    if(aggregate.method == "subsample"){
      SNPIndex = sample(pSample, p, replace = FALSE)
    } else if (aggregate.method == "bootstrap"){
      SNPIndex = sample(p, p, replace = TRUE)
      uniqueSNPIndex = SNPIndex[!duplicated(SNPIndex)]
    }
    betaHat = as.matrix(betas[SNPIndex, ], nrow = p)
    betaHat2 = as.matrix(betas[seq(pTotal + 1 - pNull, pTotal), ], nrow = pNull)

    covBetaP = betaHat %*% t(betaHat)
    covBetaQ = t(betaHat) %*% betaHat
    covBetaQ2 = t(betaHat2) %*% betaHat2

    # V
    if(V.subtract){VHat= t(svd(covBetaQ - covBetaQ2 * p / pNull)$v[,seq(1,r)])
    } else {VHat = t(svd(covBetaQ)$v[,seq(1,r)])}

    # U
    if(U.project){
      UHat = betaHat %*% t(VHat) %*% solve(VHat %*% t(VHat))
      # UHat = svd(betaHat %*% t(VHat))$u # used for metabolites data
    } else {
      # UHat = betaHat %*% t(VHat)
      UHat = svd(covBetaP)$v[,seq(1,r)] # used for metabolites data
    }

    # linear transformation
    if(rotate){
      if(is.null(parameters$rotate.method)){parameters$rotate.method = "varimax"}
      if(is.null(parameters$rotate.lambda)){parameters$rotate.lambda = 0}
      if(parameters$rotate.method == "promax"){
        if(is.null(parameters$rotate.m)){parameters$rotate.m = 4}
      }
      if(parameters$rotate.method == "varimax"){
        temp = biVarimax(U = UHat, V = VHat, lambda = parameters$rotate.lambda, normalize = TRUE)
        VHat = t(temp$loadingsV); UHat = temp$loadingsU
      } else if (parameters$rotate.method == "promax"){
        temp = biPromax(U = UHat, V = VHat, lambda = parameters$rotate.lambda, normalize = TRUE, m = parameters$rotate.m)
        VHat = t(temp$loadingsV); UHat = temp$loadingsU
      }
    }
    # if(rotate){
    #   rotate = promax(t(VHat), m = parameters$rotate.m)
    #   VHat = t(rotate$loadings)
    #   UHat = UHat %*% solve(t(rotate$rotmat))
    # }

    # truncation
    if(UV.truncate){
      UHat = myHardThresholding(UHat, zeroProp = rep(parameters$probUZero, r), normalize = FALSE)
      VHat = t(myHardThresholding(VHat, zeroProp = rep(parameters$probVZero, r), normalize = FALSE))
    }

    # post-process
    if(aggregate.method == "subsample"){
      colnames(VHat) = phenotypeNames; rownames(UHat) = SNPNames[SNPIndex]
    } else if(aggregate.method == "bootstrap"){
      UHat = UHat[!duplicated(SNPIndex),]
      colnames(VHat) = phenotypeNames; rownames(UHat) = SNPNames[uniqueSNPIndex]
    }
      VHatList[[l]] = t(VHat); UHatList[[l]] = UHat

    if(parameters$similarity.truncate){
      if(!UV.truncate){stop("similarity.truncate is only available when UV.truncate is TRUE")}
      resultDist[l,c(uniqueSNPIndex, seq(p+1, p+q)), c(uniqueSNPIndex, seq(p+1, p+q))] = rbind(UHatList[[l]] != 0, (VHatList[[l]] != 0)) %*% cbind(t(UHatList[[l]] != 0), t(VHatList[[l]]) != 0)
    }
  }

  if(parameters$similarity.truncate){
    if(aggregate.method != "bootstrap"){stop("similarity.truncate is only available for bootstrapping")}
    weightMatrixBootstrap = apply(resultDist, c(2,3), mean)
    weightMatrixBootstrap[-seq(1,p),-seq(1,p)] = (apply(resultDist[,-seq(1,p),-seq(1,p)], c(2,3), quantile, probs = max(0.20, (parameters$probVZero)^r - 1e-3)) != 0)
    weightMatrixBootstrap[seq(1,p),seq(1,p)] = (apply(resultDist[,seq(1,p),seq(1,p)], c(2,3), quantile, probs = max(0.20, 1-(1-exp(-1))^2 + sqrt((1-(1-exp(-1))^2) * (1-exp(-1))^2/parameters$n.subsample))) != 0)
    weightMatrixBootstrap[seq(1,p),-seq(1,p)] = (apply(resultDist[,seq(1,p),-seq(1,p)], c(2,3), quantile, probs = max(0.20, 1 - (1-exp(-1)) * (1-(parameters$probVZero)^r))) != 0)
    weightMatrixBootstrap[-seq(1,p),seq(1,p)] = t(weightMatrixBootstrap[seq(1,p),-seq(1,p)])
    rownames(weightMatrixBootstrap) = c(SNPNames[1:p], phenotypeNames); colnames(weightMatrixBootstrap) = c(SNPNames[1:p], phenotypeNames)
  }

  result = list(); result$VHat = VHatList; result$UHat = UHatList
  if(parameters$similarity.truncate){result$weightMatrixBootstrap = weightMatrixBootstrap}
  return(result)
}
