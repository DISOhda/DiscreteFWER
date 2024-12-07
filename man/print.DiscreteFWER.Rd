% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_fun.R
\name{print.DiscreteFWER}
\alias{print.DiscreteFWER}
\title{Printing discrete FWER results}
\usage{
\method{print}{DiscreteFWER}(x, ...)
}
\arguments{
\item{x}{object of class \code{DiscreteFWER}.}

\item{...}{further arguments to be passed to or from other methods. They
are ignored in this function.}
}
\value{
The respective input object is invisibly returned via \code{invisible(x)}.
}
\description{
Prints the results of discrete FWER analysis, stored in a \code{DiscreteFWER}
S3 class object.
}
\examples{
X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
N1 <- rep(148, 9)
N2 <- rep(132, 9)
Y1 <- N1 - X1
Y2 <- N2 - X2
df <- data.frame(X1, Y1, X2, Y2)
df

# Computation of p-values and their supports with Fisher's exact test
library(DiscreteTests)  # for Fisher's exact test
test_results <- fisher_test_pv(df)
raw_pvalues <- test_results$get_pvalues()
pCDFlist <- test_results$get_pvalue_supports()

# d-Holm with critical values; using test results object
DHolm_crit <- DHolm(test_results, critical.values = TRUE)
# print results
print(DHolm_crit)

}
