#' @examples
#' X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#' X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
#' N1 <- rep(148, 9)
#' N2 <- rep(132, 9)
#' Y1 <- N1 - X1
#' Y2 <- N2 - X2
#' df <- data.frame(X1, Y1, X2, Y2)
#' df
#'
#' # Computation of p-values and their supports with Fisher's exact test
#' library(DiscreteTests)  # for Fisher's exact test
#' test_results <- fisher_test_pv(df)
#' raw_pvalues <- test_results$get_pvalues()
#' pCDFlist <- test_results$get_pvalue_supports()
#' 
