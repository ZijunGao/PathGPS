% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{clusteringError}
\alias{clusteringError}
\title{clusteringError}
\usage{
clusteringError(U = NULL, V, UTrue = NULL, VTrue, probUZero = 0, probVZero = 0)
}
\arguments{
\item{V:}{phenotype - mediator matrix; each row represents a phenotype, each column represents a mediator. The number of mediators should be no larger than that of phenotypes.}

\item{VTrue:}{true phenotype - mediator matrix; each row represents a phenotype, each column represents a mediator. The number of mediators should be no larger than that of phenotypes.}

\item{U:}{genotype - mediator matrix estimator; each row represents a genotype, each column represents a mediator. The number of mediators should be no larger than that of genotypes. Default is NULL.}

\item{UTrue:}{true genotype - mediator matrix; each row represents a genotype, each column represents a mediator. The number of mediators should be no larger than that of genotypes. Default is NULL.}
}
\description{
compute the clustering error
}
