#' weight U statistics
#'
#' @param Y matrix of one or more response variables
#' @param G matrix of genotype
#' @param X matrix of covariate
us2 <- function(Y, G, X=NULL, acc=1e-9, ...)
{
    ## response and covariate
    Y <- as.matrix(Y)
    N <- nrow(Y)

    ## normal rank quantile standardization
    Y <- scale(apply(Y, 2, rank))

    ## covariates
    if(is.null(X) || length(X) == 0)
        X = matrix(1, N, 1L)
    else
        X = cbind(1, as.matrix(X))
    
    ## take residual of Y and G on covariates, standardize
    Y <- scale(lm(Y ~ X - 1)$residual)
    Y <- Y * sqrt((N - ncol(X)) / N)

    G <- scale(lm(G ~ X - 1)$residual)
    G <- G * sqrt((N - ncol(X)) / N)
    
    ## product kernel for genotype
    K <- tcrossprod(G) / ncol(G)

    ## product kernel for phenotype
    J <- tcrossprod(Y) / ncol(Y)

    diag(K) <- 0; # ??
    
    ## compute u score
    u <- sum(K * J)
    
    ## calculate p-value of u.
    coef <- eigen(K, symmetric=TRUE, only.values=TRUE)$values;
    pval <- imhof(u, coef)$Qq
    ## pval = davies(u, coef, acc=acc)$Qq
    pval
}


sm2 <- function(N=2e2, M=1, L=4, Q=0, goy=0, xoy=1, xog=1, eps=1, rep=1e3, ...)
{
    ret <- replicate(rep,
    {
        Y <- matrix(rnorm(N * M), N, M)
        G <- matrix(rnorm(N * L), N, L)
        X <- matrix(rnorm(N * Q), N, Q)
        
        ## confounding effect, or covariate effect only
        XOY <- X %*% matrix(rnorm(Q * M), Q, M)
        XOG <- X %*% matrix(rnorm(Q * L), Q, L)
        G <- G + XOG * xog
        Y <- Y + XOY * xoy

        ## genetic effect
        GOY <- G %*% rnorm(L * M)
        Y <- Y + GOY * goy

        us2(Y, G, X)
    })
    ret
}

ksd <- function(K)
{
    K - outer(rowMeans(K), colMeans(K), `+`) + mean(K)
}
