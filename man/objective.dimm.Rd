\name{objective.dimm}
\alias{objective.dimm}
\title{goodness-of-fit test}
\usage{
objective.dimm(object...)
}
\arguments{
\item{object}{an object of class "dimm" obtained from fitting a regression model using the dimm function}

\item{...}{additional arguments to be passed}
}
\value{
an object of class "objective.dimm". The function print is used to print the results.
}

\description{
Performs a goodness-of-fit test based on the chi-squared distribution of the quadratic objective function used by the distributed and integrated method of moments.
}
\details{
Returns the value of the quadratic form evaluated at the estimator, with corresponding degrees of freedom and p-value. objective.dimm calls the function objective.dimm.compute for actual computation.
}
