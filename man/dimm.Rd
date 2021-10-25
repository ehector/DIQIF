\name{dimm}
\alias{dimm}
\title{fit a regression model with high dimensional correlated outcome}
\usage{
dimm(formula, data, id = id, response_indicator = NULL,
  subject_indicator = NULL, family, corstr, analysis = NULL, cluster = NULL, comb_scheme=NULL...)
}
\arguments{
\item{formula}{an object of class "formula" of the type \code{Y ~ X} indicating the response ~ covariate relationship to model}

\item{data}{a data frame containing the variables in the model}

\item{id}{subject id of participants}

\item{response_indicator}{a vector of integers from 1 up to M, indicating which block the response belongs to}

\item{subject_indicator}{a vector of integers from 1 up to N, indicating which subject block the subject belongs to. Default is one block}

\item{family}{"gaussian" for continuous response, "binomial" for dichotomous response, "poisson" for count response}

\item{corstr}{correlation structure to model the covariance of the response within each block: "CS"" for compound symmetry, "AR1" for autoregressive(1)", "independence"", or CS+AR1 for a combination of compound symmetry and autoregressive(1) (only available with analysis="QIF"")}

\item{analysis}{inference method for block analysis: "PCL" for pairwise composite likelihood; "GEE" for generalized estimating equations; "QIF" for quadratic inference functions}

\item{cluster}{optional number of threads for parallelization of block analyses}

\item{comb_scheme}{optional vector of length JK specifying the combination scheme of mean parameters. For example, comb_scheme with all values set to 1 corresponds to one partition, i.e. full homogeneity}

\item{...}{additional arguments to be passed}
}
\value{
an object of class "dimm". The function summary is used to obtain and print a summary of the results. The generic accessor functions coefficients, vcov and fitted.values extract features of the object.
}
\description{
Divides the data into blocks, analyzes blocks using composite likelihood, generalized estimating equations or quadratic inference functions, and combines results using GMM
}
\details{
A formula statement has an implicit intercept. To remove the intercept, use either \code{Y ~ X-1} or \code{X ~ 0+X}. dimm calls functions dimm.compute.mean or dimm.compute.mean.parallel for actual computation.
}
