# helper of metabolites

#' myDistance
#'
#' compute the distance between the col space of X1 and that of X2
#' @param X covariate matrix to predict conditional densities at; each row represents an observation vector.
#'
#' @export
myDistance = function(X1,X2){
  if(is.null(dim(X1))){X1 = matrix(X1, ncol = 1); X2 = matrix(X2, ncol = 1)}
  if(is.null(dim(X2))){X1 = matrix(X1, ncol = 1); X2 = matrix(X2, ncol = 1)}
  if(rankMatrix(t(X1) %*% X1)[1] < dim(X1)[2]){dist = 0; return (dist)}
  if(rankMatrix(t(X2) %*% X2)[1] < dim(X2)[2]){dist = 0; return (dist)}
  P1 = X1 %*% solve(t(X1) %*% X1) %*% t(X1)
  P2 = X2 %*% solve(t(X2) %*% X2) %*% t(X2)
  dist = 1 - max(abs(eigen(P1-P2)$values))
  return (dist)
}

# hard thresholding for sparse estimation
# input:
  # V: p by q vector, we truncate each column and preserve proportions of zeros in each column at zeroProp
  # zeroProp:proprotion of zeros for each column; if only a scalar is input, control the proportion of zeros in V

myHardThresholding = function(V, zeroProp){
  if(is.null(V)){V = matrix(V, ncol = 1)}
  p = dim(V)[1]; q = dim(V)[2]
  if(length(zeroProp) == 1){
    # correct: at least preserve the largest...
    index1 = apply(abs(V),2,which.max)
    index = order(abs(V), decreasing = TRUE)[1:floor(zeroProp * p * q)]
    V[-index] = 0; V = matrix(V, ncol = q)
    VNorm = diag(t(V) %*% V)
    V[cbind(index1[VNorm == 0], which(VNorm == 0))] = 1; VNorm[which(VNorm == 0)] = 1
    V = V %*% diag(1/sqrt(VNorm))
    return(V)
  }
  for(j in 1:q){
    index = order(abs(V[,j]), decreasing = TRUE)[1:floor(zeroProp[j] * p)]
    V[-index, j] = 0
    V[,j] = V[,j]/sqrt(sum((V[,j])^2))
  }
  return(V)
}


# transfrorm EVs to a weighted undirected graph
EVToGraph = function(V, varProp = NULL, k = NULL, eigenvalues = NULL, returnFull = FALSE){
  edgeList = list()
  if(class(V) != "list"){
    r = min(dim(V))
    edgeList = EVToGraphHelper(V = V, edgeList = edgeList,
                               varProp = varProp, k = k, eigenvalues = eigenvalues)
  } else {
    r = min(dim(V[[1]]))
    for(i in 1:length(V)){
      edgeList = EVToGraphHelper(V = V[[i]], edgeList = edgeList,
                                 varProp = varProp, k = k, eigenvalues = eigenvalues)
    }
  }
  # construct graphs
  if(length(edgeList) == 0){print("no edge"); return(NULL)}
  weight = table(unlist(edgeList))/r; if(class(V) == "list"){weight = weight/length(V)}
  edge = unlist(strsplit(names(weight),split = ";")); g = igraph::graph(edge) # isolates not needed
  names(weight) = NULL; igraph::E(g)$weight = weight
  g = igraph::as.undirected(g, mode= "collapse",edge.attr.comb=list(weight="sum", "ignore"))
  if(!returnFull){return(g)
    } else {result = list(); result$g = g; result$edge = edge; result$weight = weight; return(result)}
}

EVToGraphHelper = function(V, edgeList, varProp = NULL, k = NULL, eigenvalues = NULL){
  # preprocess
  if(is.null(dim(V))){V = matrix(V, ncol = 1)}
  if(dim(V)[1] <= dim(V)[2]){V = t(V)}
  q = dim(V)[1]; r = dim(V)[2]
  if(is.null(rownames(V))){stop("please input names of V")} else {nodeNames = rownames(V)}
  edgeCounter = length(edgeList)
  if(is.null(eigenvalues)){
    # determine dominating phenotypes
    if(!is.null(k)){
      for(i in 1:r){
        index = order((V[,i])^2, decreasing = TRUE)[1:k]
        for(j in 2:k){
          for(l in 1:(j-1)){
            edgeCounter = edgeCounter + 1
            edgeList[[edgeCounter]] = paste(nodeNames[index[j]], nodeNames[index[l]], sep = ";")
          }
        }
      }
    } else if(!is.null(varProp)){
      for(i in 1:r){
        maxIndex = suppressWarnings(max(which(cumsum(sort((V[,i])^2, decreasing = TRUE)) <= varProp)))
        if(maxIndex <= 1){next}
        index = order((V[,i])^2, decreasing = TRUE)[1:maxIndex]
        for(j in 2:maxIndex){
          for(l in 1:(j-1)){
            edgeCounter = edgeCounter + 1
            edgeList[[edgeCounter]] = paste(nodeNames[index[j]], nodeNames[index[l]], sep = ";")
          }
        }
      }
    }
  }
  return (edgeList)
}

weightToGraph = function(weightMatrix, mode = "undirected"){
  diag(weightMatrix) = 0
  g = igraph::graph_from_adjacency_matrix(weightMatrix, mode = mode, weighted = TRUE)
  g = igraph::as.undirected(g, mode= "collapse",edge.attr.comb=list(weight="sum", "ignore"))
  return(g)
}

# transfrorm EVs to a full weighted undirected graph with three types of edges
EVToGraphFull = function(U, V, varPropU = NULL, kU = NULL, varPropV = NULL, kV = NULL, betweenWeight = 1, phenotypeWeight = 1, genotypeWeight = 1, eigenvalues = NULL, returnFull = FALSE){
  edgeList = list(); edgeList$edgeListUU = list(); edgeList$edgeListVV = list(); edgeList$edgeListUV = list()
  if(class(V) != "list"){
    r = min(dim(V))
    edgeList = EVToGraphFullHelper(U = U, V = V, edgeList = edgeList,
                               varPropU = varPropU, kU = kU, varPropV = varPropV,
                               kV = kV, eigenvalues = eigenvalues)
  } else {
    r = min(dim(V[[1]]))
    for(i in 1:length(V)){
      edgeList = EVToGraphFullHelper(U = U[[i]], V = V[[i]], edgeList = edgeList,
                                    varPropU = varPropU, kU = kU, varPropV =
                                    varPropV, kV = kV, eigenvalues = eigenvalues)
    }
  }
  # construct graphs
  if((length(edgeList$edgeListUU) + length(edgeList$edgeListVV) + length(edgeList$edgeListUV)) == 0){print("no edge"); return(NULL)}
  weight = c(table(unlist(edgeList$edgeListUU)) * genotypeWeight, table(unlist(edgeList$edgeListVV)) * phenotypeWeight, table(unlist(edgeList$edgeListUV)) * betweenWeight)
  indexZeroWeight = which(weight == 0)
  if(length(indexZeroWeight) != 0){weight = weight[-indexZeroWeight]}
  weight = weight/r; if(class(V) == "list"){weight = weight/length(V)}
  edge = unlist(strsplit(names(weight),split = ";")); g = igraph::graph(edge) # isolates not needed
  names(weight) = NULL; igraph::E(g)$weight = weight
  g = igraph::as.undirected(g, mode= "collapse",edge.attr.comb=list(weight="sum", "ignore"))
  if(!returnFull){return(g)
  } else {
    result = list(); result$g = g; result$edge = edge; result$weight = weight
    return(result)
  }
}

EVToGraphFullHelper = function(U, V, edgeList, varPropU = NULL, kU = NULL, varPropV = NULL, kV = NULL, eigenvalues = NULL){
  # preprocess
  if(is.null(dim(V))){V = matrix(V, ncol = 1)}
  if(is.null(dim(V))){U = matrix(U, ncol = 1)}
  if(dim(V)[1] <= dim(V)[2]){V = t(V)}
  if(dim(U)[1] <= dim(U)[2]){U = t(U)}
  q = dim(V)[1]; r = dim(V)[2]; p = dim(U)[1]
  if(is.null(rownames(V))){stop("please input names of V")} else {nodeNamesV = rownames(V)}
  if(is.null(rownames(U))){stop("please input names of U")} else {nodeNamesU = rownames(U)}
  edgeCounterVV = length(edgeList$edgeListVV)
  edgeCounterUV = length(edgeList$edgeListUV)
  edgeCounterUU = length(edgeList$edgeListUU)
  if(is.null(eigenvalues)){
    # determine dominating phenotypes
    # VV
    if(!is.null(kV)){
      for(i in 1:r){
        indexV = order((V[,i])^2, decreasing = TRUE)[1:kV]
        for(j in 2:kV){
          for(l in 1:(j-1)){
            edgeCounterVV = edgeCounterVV + 1
            edgeList$edgeListVV[[edgeCounterVV]] = paste(nodeNamesV[indexV[j]], nodeNamesV[indexV[l]], sep = ";")
          }
        }
      }
    } else if(!is.null(varPropV)){
      for(i in 1:r){
        maxIndexV = suppressWarnings(max(which(cumsum(sort((V[,i])^2, decreasing = TRUE)) <= varPropV)))
        if(maxIndexV <= 1){next} else {
          indexV = order((V[,i])^2, decreasing = TRUE)[1:maxIndexV]
          for(j in 2:maxIndexV){
            for(l in 1:(j-1)){
              edgeCounterVV = edgeCounterVV + 1
              edgeList$edgeListVV[[edgeCounterVV]] = paste(nodeNamesV[indexV[j]], nodeNamesV[indexV[l]], sep = ";")
            }
          }
        }
      }
    }
    # UU
    if(!is.null(kU)){
      for(i in 1:r){
        indexU = order((U[,i])^2, decreasing = TRUE)[1:kU]
        for(j in 2:kU){
          for(l in 1:(j-1)){
            edgeCounterUU = edgeCounterUU + 1
            edgeList$edgeListUU[[edgeCounterUU]] = paste(nodeNamesU[indexU[j]], nodeNamesU[indexU[l]], sep = ";")
          }
        }
      }
    } else if(!is.null(varPropU)){
      for(i in 1:r){
        maxIndexU = suppressWarnings(max(which(cumsum(sort((U[,i])^2, decreasing = TRUE)) <= varPropU)))
        if(maxIndexU <= 1){next}
        indexU = order((U[,i])^2, decreasing = TRUE)[1:maxIndexU]
        for(j in 2:maxIndexU){
          for(l in 1:(j-1)){
            edgeCounterUU = edgeCounterUU + 1
            edgeList$edgeListUU[[edgeCounterUU]] = paste(nodeNamesU[indexU[j]], nodeNamesU[indexU[l]], sep = ";")
          }
        }
      }
    }
    # UV
    if((is.null(kU) + is.null(kV)) == 0){
      for(i in 1:r){
        indexU = order((U[,i])^2, decreasing = TRUE)[1:kU]
        indexV = order((V[,i])^2, decreasing = TRUE)[1:kV]
        for(j in 1:kU){
          for(l in 1:kV){
            edgeCounterUV = edgeCounterUV + 1
            edgeList$edgeListUV[[edgeCounterUV]] = paste(nodeNamesU[indexU[j]], nodeNamesV[indexV[l]], sep = ";")
          }
        }
      }
    } else if((is.null(kV) + is.null(varPropU)) == 0){
      for(i in 1:r){
        maxIndexU = suppressWarnings(max(which(cumsum(sort((U[,i])^2, decreasing = TRUE)) <= varPropU)))
        indexV = order((V[,i])^2, decreasing = TRUE)[1:kV]
        if(maxIndexU <= 0){next}
        indexU = order((U[,i])^2, decreasing = TRUE)[1:maxIndexU]
        for(j in 1:maxIndexU){
          for(l in 1:kV){
            edgeCounterUV = edgeCounterUV + 1
            edgeList$edgeListUV[[edgeCounterUV]] = paste(nodeNamesU[indexU[j]], nodeNamesV[indexV[l]], sep = ";")
          }
        }
      }
    } else if ((is.null(varPropV) + is.null(kU)) == 0){
      for(i in 1:r){
        indexU = order((U[,i])^2, decreasing = TRUE)[1:kU]
        maxIndexV = suppressWarnings(max(which(cumsum(sort((V[,i])^2, decreasing = TRUE)) <= varPropV)))
        if(maxIndexV <= 0){next}
        indexV = order((V[,i])^2, decreasing = TRUE)[1:maxIndexV]
        for(j in 1:kU){
          for(l in 1:maxIndexV){
            edgeCounterUV = edgeCounterUV + 1
            edgeList$edgeListUV[[edgeCounterUV]] = paste(nodeNamesU[indexU[j]], nodeNamesV[indexV[l]], sep = ";")
          }
        }
      }
    } else if ((is.null(varPropV) + is.null(varPropU)) == 0){
      for(i in 1:r){
        maxIndexU = suppressWarnings(max(which(cumsum(sort((U[,i])^2, decreasing = TRUE)) <= varPropU)))
        maxIndexV = suppressWarnings(max(which(cumsum(sort((V[,i])^2, decreasing = TRUE)) <= varPropV)))
        if(maxIndexU <= 0){next}
        indexU = order((U[,i])^2, decreasing = TRUE)[1:maxIndexU]
        indexV = order((V[,i])^2, decreasing = TRUE)[1:maxIndexV]
        for(j in 1:maxIndexU){
          for(l in 1:maxIndexV){
            edgeCounterUV = edgeCounterUV + 1
            edgeList$edgeListUV[[edgeCounterUV]] = paste(nodeNamesU[indexU[j]], nodeNamesV[indexV[l]], sep = ";")
          }
        }
      }
    }
  }
  return (edgeList)
}


# bipartite graph
EVToGraphBipartite = function(V, varProp = NULL, k = NULL, eigenvalues = NULL, returnFull = FALSE){
  if(class(V) != "list"){
    edgeList = EVToGraphHelper(V = V, edgeList = list(),
                               varProp = varProp, k = k, eigenvalues = eigenvalues)
  } else {
    edgeList = list()
    for(i in 1:length(V)){
      edgeList = EVToGraphHelper(V = V[[i]], edgeList = edgeList,
                                 varProp = varProp, k = k, eigenvalues = eigenvalues)
    }
  }
  # construct graphs
  if(length(edgeList) == 0){print("no edge"); return(NULL)}
  weight = table(unlist(edgeList))/r; if(class(V) == "list"){weight = weight/length(V)}
  edge = unlist(strsplit(names(weight),split = ";")); g = igraph::graph(edge) # isolates not needed
  names(weight) = NULL; igraph::E(g)$weight = weight
  g = igraph::as.undirected(g, mode= "collapse",edge.attr.comb=list(weight="sum", "ignore"))
  if(!returnFull){return(g)
  } else {result = list(); result$g = g; result$edge = edge; result$weight = weight; return(result)}
}

EVToGraphBipartiteHelper = function(V, edgeList, varProp = NULL, k = NULL, eigenvalues = NULL){
  # preprocess
  if(is.null(dim(V))){V = matrix(V, ncol = 1)}
  if(dim(V)[1] <= dim(V)[2]){V = t(V)}
  q = dim(V)[1]; r = dim(V)[2]
  if(is.null(rownames(V))){stop("please input names of V")} else {nodeNames = rownames(V)}
  edgeCounter = length(edgeList)
  if(is.null(eigenvalues)){
    # determine dominating phenotypes
    if(!is.null(k)){
      for(i in 1:r){
        index = order((V[,i])^2, decreasing = TRUE)[1:k]
        for(j in 2:k){
          for(l in 1:(j-1)){
            edgeCounter = edgeCounter + 1
            edgeList[[edgeCounter]] = paste(nodeNames[index[j]], nodeNames[index[l]], sep = ";")
          }
        }
      }
    } else if(!is.null(varProp)){
      for(i in 1:r){
        maxIndex = suppressWarnings(max(which(cumsum(sort((V[,i])^2, decreasing = TRUE)) <= varProp)))
        if(maxIndex <= 1){return (edgeList)}
        index = order((V[,i])^2, decreasing = TRUE)[1:maxIndex]
        for(j in 2:maxIndex){
          for(l in 1:(j-1)){
            edgeCounter = edgeCounter + 1
            edgeList[[edgeCounter]] = paste(nodeNames[index[j]], nodeNames[index[l]], sep = ";")
          }
        }
      }
    }
  }
  return (edgeList)
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

