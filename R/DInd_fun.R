#' @name DInd
#' 
#' @title
#' Discrete FWER Procedure for Independent Tests
#' 
#' @description 
#' `DInd()` is a wrapper function of [`discrete_FWER()`] for computing discrete
#' MTP adaptations for discrete tests. It simply passes its arguments to
#' [`discrete_FWER()`] with fixed `independence = TRUE`.
#' 
#' @templateVar test_results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar stepdown TRUE
#' @templateVar critical_values TRUE
#' @templateVar select_threshold TRUE
#' @templateVar pCDFlist_indices TRUE
#' @templateVar triple_dots TRUE
#' @template param
#' 
#' @template details_crit
#' 
#' @template return
#' 
#' @seealso
#' [`discrete_FWER()`], [`DBonf()`], [`DHolm()`]
#' 
#' @references
#' Döhler, S. (2010). Validation of credit default probabilities using
#'   multiple-testing procedures. *Journal of Risk Model Validation*, *4*(4),
#'   59-92. \doi{10.21314/JRMV.2010.062}
#' 
#' @template example
#' @examples
#' # d-Ind without critical values; using extracted p-values and supports
#' DInd_fast <- DInd(raw_pvalues, pCDFlist)
#' summary(DInd_fast)
#' 
#' # d-Ind with critical values; using test results object
#' DInd_crit <- DInd(test_results, critical_values = TRUE)
#' summary(DInd_crit)
#'  
#' # d-Ind (step-down) without critical values; using test results object
#' DInd_sd_fast <- DInd(test_results, stepdown = TRUE)
#' summary(DInd_sd_fast)
#' 
#' # d-Ind (step-down) with critical values; using extracted p-values and supports
#' DInd_sd_crit <- DInd(raw_pvalues, pCDFlist, stepdown = TRUE, critical_values = TRUE)
#' summary(DInd_sd_crit)
#' 
#' @export
DInd <- function(test_results, ...) UseMethod("DInd")

#' @rdname DInd
#' @export
DInd.default <- function(
    test_results,
    pCDFlist,
    alpha            = 0.05,
    stepdown         = FALSE,
    critical_values  = FALSE,
    select_threshold = 1,
    pCDFlist_indices = NULL,
    ...
) {
  out <- discrete_FWER.default(
    test_results     = test_results,
    pCDFlist         = pCDFlist,
    alpha            = alpha,
    independence     = TRUE,
    stepdown         = stepdown,
    critical_values  = critical_values,
    select_threshold = select_threshold,
    pCDFlist_indices = pCDFlist_indices,
    ...
  )
  
  out$Data$Data_name <- paste(
    deparse(substitute(test_results)),
    "and",
    deparse(substitute(pCDFlist))
  )
  
  return(out)
}

#' @rdname DInd
#' @export
DInd.DiscreteTestResults <- function(
    test_results,
    alpha            = 0.05,
    stepdown         = FALSE,
    critical_values  = FALSE,
    select_threshold = 1,
    ...
) {
  out <- discrete_FWER.DiscreteTestResults(
    test_results     = test_results,
    alpha            = alpha,
    independence     = TRUE,
    stepdown         = stepdown,
    critical_values  = critical_values,
    select_threshold = select_threshold,
    ...
  )
  
  out$Data$Data_name <- deparse(substitute(test_results))
  
  return(out)
}
