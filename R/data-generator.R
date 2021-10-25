posdef.matrix <- function(n, seed){
  ## Function written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708
  set.seed(seed)
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  O <- Q %*% diag(d / abs(d))
  Z <- t(O) %*% diag(runif(n, 0, 10)) %*% O
  return(Z)
}

#' @name dataset.normal.X2
#' @title simulate data
#' @description Simulate data for illustration of DIQIF
#' @param family "gaussian" for continuous response
#' @param N sample size
#' @param M dimension of response
#' @param theta vector of true parameter values 
#' @param response_indicator a vector of integers from 1 up to M, indicating which block the response belongs to
#' @param seed set random seed for data generation
#' @keywords data simulator
#' @details Simulates an artificial dataset with two covariates and response with block AR-1 covariance structure. 
#' Covariate effects (theta) are assumed homogeneous across outcome groups for simplicity.
#' @return returns a simulated dataset with subject id, response, and two covariates X1 and X2
#' @export
#' 
dataset.normal.X2 <- function(family, N, M, theta, response_indicator, seed){
  set.seed(seed)
  sd <- 2.0
  r <- 0.5
  A <- sd^2 * r^abs(outer(1:max(table(response_indicator)), 1:max(table(response_indicator)) , "-"))
  S <- posdef.matrix(length(unique(response_indicator)), seed=seed)
  
  if (length(unique(table(response_indicator)))==1){
    Sigma <- kronecker(S,A)
  } else{
    Sigma <- kronecker(S,A)
    for(i in 1:length(unique(response_indicator))){
      if(table(response_indicator)[i] < max(table(response_indicator))){
        start <- sum(table(response_indicator)[1 : i])+1
        end <- start + max(table(response_indicator))-table(response_indicator)[i]-1
        Sigma <- Sigma[-c(start:end), -c(start:end)]
      }
    }
  }
  
  ## Generate N replicates of an M-dimensional Normal vector with mean theta[1]+theta[2]*x1+theta[3]*x2
  ## (X_1,X_2) ~ N_M, Y ~ N_M(theta[1]+theta[2]*x1+theta[3]*x2, Sigma)
  Sigma_list <- as.matrix(Matrix::bdiag(lapply(unique(subject_indicator), function(x) { 
    set.seed(seed+1000*x)
    Sigma*runif(1,0,1) 
  })))
  set.seed(seed)
  epsilon <- MASS::mvrnorm(N/length(unique(subject_indicator)), rep(0,M*length(unique(subject_indicator))), Sigma_list, tol = 1e-6)
  covariates <- MASS::mvrnorm(N, c(rep(0,M), rep(1,M)), as.matrix(Matrix::bdiag(posdef.matrix(M, seed=seed), posdef.matrix(M, seed=seed*2000))))
  X_1 <- covariates[,1:M]
  X_2 <- covariates[,(M+1):(2*M)]
  new_eps <- do.call(rbind, lapply(unique(subject_indicator), function(x){ epsilon[,((x-1)*M+1):(((x-1)*M+M))] }))
  Y <- t(new_eps + rep(theta[1], M)+theta[2]*X_1+theta[3]*X_2)
  
  data_short <- cbind(id=seq(1,N), X1=X_1, X2=X_2)
  colnames(data_short)[2:(M+1)] <- paste("X1", seq(1, M, 1), sep="_")
  colnames(data_short)[(M+2):(2*M+1)] <- paste("X2", seq(1, M, 1), sep="_")
  data <- reshape(as.data.frame(data_short), direction="long", varying=c(2:(2*M+1)), sep="_")
  data <- data[order(data$id),-which(names(data)=="time")]
  data <- cbind(data, response=c(Y))
  return(data)
}
