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
#' @param m
#' @param V.subtract
#' @param V.rotate
#' @param V.rotate.m
#' @param U.project
#' @param randomSeed
#'
#' @return VHat, UHat
#'
#' @export
columnEstimator = function(betas = "metabolites", p, pNull, pSample = NULL, r, m = 100, V.subtract = TRUE, V.rotate = TRUE, V.rotate.m = 1, U.project = TRUE, randomSeed = 318){
  # preprocess
  if(is.character(betas)){
    if(betas == "metabolites"){
      betas = read.table(file = "./data/metabolitesBeta.txt", header = TRUE)
      suppressWarnings(RNGkind(sample.kind = "Rounding"))
    }
  } else if(is.null(dim(betas))){stop("please input a matrix of summary statistics for multiple genotypes and multiple phenotypes")}
  pTotal = dim(betas)[1]; q = dim(betas)[2]
  if(is.null(pSample)){pSample = 2 * p}

  if(!is.null(colnames(betas))){phenotypeNames = colnames(betas)
  } else {phenotypeNames = paste("phenotype", seq(1,q), sep = "")}
  if(!is.null(rownames(betas))){SNPNames = rownames(betas)
  } else {SNPNames = paste("genotype", seq(1,pTotal), sep = "")}

  set.seed(randomSeed)

  # estimation
  VHatList = list(); UHatList = list()
  for(l in 1:m){
    SNPIndex = sample(pSample, p, replace = FALSE)
    betaHat = as.matrix(betas[SNPIndex, ], nrow = p)
    betaHat2 = as.matrix(betas[seq(pTotal + 1 - pNull, pTotal), ], nrow = pNull)

    covBetaP = betaHat %*% t(betaHat)
    covBetaQ = t(betaHat) %*% betaHat
    covBetaQ2 = t(betaHat2) %*% betaHat2

    # V
    if(V.subtract){VHat= t(svd(covBetaQ - covBetaQ2 * p / pNull)$v[,seq(1,r)])
    } else {VHat = t(svd(covBetaQ)$v[,seq(1,r)])}

    # U
    if(U.project){UHat = svd(betaHat %*% t(VHat))$u
    } else {UHat = svd(covBetaP)$v[,seq(1,r)]}

    # linear transformation
    if(V.rotate){
      rotate = promax(t(VHat), m = V.rotate.m)
      VHat = t(rotate$loadings)
      UHat = UHat %*% solve(t(rotate$rotmat))
    }

    # post-process
    colnames(VHat) = phenotypeNames; rownames(UHat) = SNPNames[SNPIndex]
    VHatList[[l]] = VHat; UHatList[[l]] = UHat
  }

  result = list(); result$VHat = VHatList; result$UHat = UHatList
  return(result)
}
