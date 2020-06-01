library(MASS)
library(tidyverse)

var_and_covar <- function(data){
    ## --- Compute variance-covariance matrix from input data
    ## Input: list of vector with the same lengths.
    ## Output: dxd Matrix, where d is the length of vectors.
    diad <- function(x) x %*% t(x)
    Reduce("+", lapply(data,diad)) / length(data)
}

distance <- function(x,y) {
    ## --- Compute 2-norm of x-y.
    ## Input: Two matrices with the same dimension.
    ## Output: Scalar.
    max(abs(eigen(x - y)$values))
}

simulate <- function(n,d,sigma){
    ## --- Sample n-dimensional normal distribution with 0-mean and sigma-variance.
    ## Input: number of the sample, the dimension of the distribution, and the variance-coveriance matrix.
    ## Output: List of vector.
    data <- mvrnorm(n,mu=rep(0,d), Sigma = sigma)
    split(data,row(data))
}

perf <- function(n,d,sigma){
    ## --- Performance check.
    data <- simulate(n,d,sigma)
    vars <- var_and_covar(data)
    distance(sigma, vars)
}

random_sigma <- function(d,p){
    ## --- Create a covariance matrix.
    ## Input: the dimension of the distribution and the probability that
    ## an element of Sigma would be 1, i.e., we compute Sigma s.t.,
    ## Sigma_{i,j} = 1 with probability p(i <= j).
    sigma <- matrix(sample(x=c(0,1), size = d * d, prob=c(p,1-p), replace=TRUE),
                    nrow = d,
                    ncol = d)
    print(sigma)
    sigma[lower.tri(sigma)] <- 0
    sigma <- sigma + t(sigma)
    diag(sigma) <- diag(sigma)/ 2
    sigma
}


## Simulate digonal case.
dimensions <- seq(from = 100, to = 20000, by = 500)
sample_num <- rep(10000,length(dimensions))
sigmas <- mapply(function(sn,d) perf(sn,d,diag(d)), sample_num, dimensions)
