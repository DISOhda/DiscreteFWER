#' @return
#' A \code{DiscreteFWER} S3 class object whose elements are:
#' \item{Rejected}{rejected raw \eqn{p}-values.}
#' \item{Indices}{indices of rejected hypotheses.}
#' \item{Num_rejected}{number of rejections.}
#' \item{Adjusted}{adjusted \eqn{p}-values.}
#' \item{Critical_constants}{critical values (only exists if computations where performed with `critical_values = TRUE`).}
#' \item{Data}{list with input data.}
#' \item{Data$Method}{character string describing the performed algorithm, e.g. 'Discrete Bonferroni procedure'.}
#' \item{Data$Raw_pvalues}{observed \eqn{p}-values.}
#' \item{Data$pCDFlist}{list of the \eqn{p}-value supports.}
#' \item{Data$FWER_level}{FWER level `alpha`.}
#' \item{Data$Independence}{boolean indicating whether the \eqn{p}-values were considered as independent.}
#' \item{Data$Stepdown}{boolean indicating whether a step-down or single-step procedure was performed.}
#' \item{Data$Data_name}{the respective variable names of the input data.}
#' \item{Select}{list with data related to \eqn{p}-value selection; only exists if `select_threshold < 1`.}
#' \item{Select$Threshold}{\eqn{p}-value selection threshold (`select_threshold`).}
#' \item{Select$Effective_Thresholds}{results of each \eqn{p}-value CDF evaluated at the selection threshold.}
#' \item{Select$Pvalues}{selected \eqn{p}-values that are \eqn{\leq} selection threshold.}
#' \item{Select$Indices}{indices of \eqn{p}-values \eqn{\leq} selection threshold.}
#' \item{Select$Scaled}{scaled selected \eqn{p}-values.}
#' \item{Select$Number}{number of selected \eqn{p}-values \eqn{\leq} selection threshold.}
