\name{dataset.normal.X2}
\alias{dataset.normal.X2}
\title{simulate data}
\usage{
dataset.normal.X2(family, N, M, theta, response_indicator, seed)
}
\arguments{
\item{family}{"gaussian" for continuous response}

\item{N}{sample size}

\item{M}{dimension of response}

\item{theta}{vector of true parameter values}

\item{response_indicator}{a vector of integers from 1 up to M, indicating which block the response belongs to}

\item{seed}{set random seed for data generation}
}
\value{
returns a simulated dataset with subject id, response, and two covariates X1 and X2
}
\description{
Simulate data for illustration of DIQIF
}
\details{
Simulates an artificial dataset with two covariates and response with block AR-1 covariance structure. Covariate effects (theta) are assumed homogeneous across outcome groups for simplicity.
}
\keyword{data}
\keyword{simulator}
