library(MASS)
library(tidyverse)
var_and_covar <- function(data){
    ## --- Compute variance-covariance matrix from input data
    ## Input: list of vector with the same lengths.
    ## Output: dxd Matrix, where d is the length of vectors.
    diad <- function(x) x %*% t(x)
    Reduce("+", lapply(data,diad)) / length(data)
}

filter_small <- function(x, lambda){
    apply(X=x,MARGIN=c(1,2),function(x)ifelse(x < lambda, 0, x))
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

perf_filter <- function(n,d,sigma){
    ## --- Performance check.
    data <- simulate(n,d,sigma)
    lambda <- sqrt(2*log(d)/n)
    vars <- filter_small(var_and_covar(data), lambda)
    distance(sigma, vars)
}


random_sigma <- function(d,p){
    ## --- Create a covariance matrix.
    ## Input: the dimension of the distribution and the probability that
    ## an element of Sigma would be 1, i.e., we compute Sigma s.t.,
    ## Sigma_{i,j} = 1 with probability p(i <= j).
    p <- sqrt(p)
    sigma <- matrix(sample(x=c(1,0), size = d * d, prob=c(p,1-p), replace=TRUE),
                    nrow = d,
                    ncol = d)
    diag(sigma) <- rep(1, d)
    sigma * t(sigma)
}


## Simulate digonal case.
start_dim <- 10
end_dim <- 200
by <- 50
replicate <- 50
sample_number <- 1000

start_sample_number <- 100
end_sample_number <- 1500
frac <- 0.2
sample_span_by <- 300

run_simulation <- function(output_name, estimate){
    dimensions <- rep(seq(from = start_dim, to = end_dim, by = by), each = replicate)
    sample_num <- rep(sample_number,length(dimensions))
    result <- tibble(
        dimension = dimensions,
        sample_num = sample_num,
        distance = mapply(estimate, sample_num, dimensions)
    )
    g <- result %>% ggplot() + geom_violin(aes(x = factor(dimension), y = distance))
    ggsave(filename = paste0(output_name, ".png"), g)
    temp <- result

    sample_num <- rep(seq(from = start_sample_number, to = end_sample_number, by = sample_span_by), each = replicate)
    dimensions <- sample_num * 0.2
    result <- tibble(
        dimension = dimensions, sample_num = sample_num,
        distance = mapply(estimate, sample_num, dimensions)
    )
    g <- result %>% ggplot() +
        geom_violin(aes(x = factor(sample_num), y = distance))
    ggsave(filename = paste0(output_name,"_sample.png"), g)
    list(first = temp, second = result)
}

result <- run_simulation("naive", function(sn, d)perf(sn,d,diag(d)))
result <- run_simulation("filter", function(sn, d)perf_filter(sn, d, diag(d)))
prob <- 0.1

result <- run_simulation("naive_additional", function(sn, d) perf(sn, d, random_sigma(d,prob)))
result <- run_simulation("filter_additional", function(sn, d) perf_filter(sn, d, random_sigma(d,prob)))
