% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bootstrap.R
\name{lir_bootstrap}
\alias{lir_bootstrap}
\title{Nonparametric bootstrapping for LIR models
observed_individual}
\usage{
lir_bootstrap(X, tp, seed = NULL)
}
\arguments{
\item{X}{A list or matrix containing the identities of individuals identified during each sampling period}

\item{tp}{A set of observed time.}

\item{seed}{Random seed}
}
\value{
The bootstrap samples.
}
\description{
Nonparametric bootstrapping for LIR models
observed_individual
}
\details{
See Efron and Tibshirani (1993) for details on this function.
}
\examples{

# Example
# generate the original data
tp <- c(1:5, 51:55, 101:105, 501:505, 601:605)
N <- 100; n <- 10
# X is a list
X <- list()
for (i in tp){
X[[i]] <- sample(1:N, n)
}
# bootstrap sample data
new_sample <- lir_bootstrap(X, tp)
# X is a matrix
matrix_data <- list_to_matrix(X, tp)
new_sample <- lir_bootstrap(X, tp)

}
\references{
Efron, B. and Tibshirani, R. (1986). The bootstrap method for
standard errors, confidence intervals, and other measures of statistical accuracy.
Statistical Science, Vol 1., No. 1, pp 1-35.

Efron, B. and Tibshirani, R. (1993) An Introduction to the Bootstrap.
Chapman and Hall, New York, London.
}
