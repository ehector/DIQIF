\name{DIQIF-package}
\alias{DIQIF-package}
\alias{DIQIF}
\docType{package}
\title{\packageTitle{DIQIF}}
\description{\packageDescription{DIQIF}}
\details{
  THIS SOURCE CODE IS SUPPLIED "AS IS"" WITHOUT WARRANTY OF ANY KIND, AND ITS AUTHOR DISCLAIMS ANY AND ALL WARRANTIES, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, AND ANY WARRANTIES OR NON INFRINGEMENT. THE USER ASSUMES ALL LIABILITY AND RESPONSIBILITY FOR USE OF THIS SOURCE CODE, AND THE AUTHOR WILL NOT BE LIABLE FOR DAMAGES OF ANY KIND RESULTING FROM ITS USE. Without limiting the generality of the foregoing, the author does not warrant that the Source Code will be error-free, will operate without interruption, or will meet the needs of the user.

  The DESCRIPTION file: \packageDESCRIPTION{DIQIF}
  \packageIndices{DIQIF}
}
\author{
  \packageAuthor{DIQIF}
  
  Maintainer: \packageMaintainer{DIQIF}
}
\references{
  The posdef.matrix function was written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.
  GEE implementation is through the R package geepack: Højsgaard, S., Halekoh, U., Yan, J. (2006). The R package geepack for generalized estimating equations. Journal of Statistical Software, 15(2):1–11.
  Emily C. Hector and Peter X.-K. Song. A distributed and integrated method of moments for high-dimensional correlated data analysis. Journal of the American Statistical Association, pages 1–14, 2020. doi: 10.1080/01621459.2020.1736082.
  
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
  \code{\link[<geepack>:<geepack>-geepack]{<geepack>}}
}
\examples{
  ########### Generate Data ########### 
N <- 1000                   # sample size
M <- 200                    # response dimension
family <- "gaussian"        # response is Gaussian distributed
method <- "exact"           # find one-step estimator
analysis <- "QIF"           # using QIF analysis method
theta <- c(0.3, 0.6, 0.8)   # the true parameter value

# The response indicator partitions the response into blocks
# The subject indicator partitions subjects into blocks
response_indicator <- c(rep(1,45), rep(2,42), rep(3,50), rep(4,34), rep(5,29))
subject_indicator <- c(rep(1,N/2), rep(2,N/2))

# Generate a dataset with two covariates, response of 
# dimension M with block AR-1 covariance structure
data <- dataset.normal.X2(family, N, M, theta, response_indicator, seed=500)

########### Analysis with DIQIF ###########

# Model block covariance with independence structure
results.ind <- dimm(response~X1+X2, data, id="id", response_indicator, subject_indicator, 
                    family, corstr="independence", analysis, cluster=5)

# Model block covariance with compound symmetry structure
results.CS <- dimm(response~X1+X2, data, id="id", response_indicator, subject_indicator, 
                   family, corstr="CS", analysis, cluster=5)

# Model block covariance with AR-1 structure                             
results.AR1 <- dimm(response~X1+X2, data, id="id", response_indicator, subject_indicator, 
                    family, corstr="AR1", analysis, cluster=5)
                                        
## Goodness-of-fit test
objective.dimm(results.AR1)

# Fit model combining first three response groups and last two response groups across subject groups                          
results.AR1.sel <- dimm(response~X1+X2, data, id="id", response_indicator, subject_indicator,
                        family, corstr="AR1", analysis, cluster=1, comb_scheme=c(1,1,1,2,2,1,1,1,2,2))
}
