#include "helper.h"

//' @name kernel
//' 
//' @keywords internal
//' 
//' @title
//' Kernel Functions
//' 
//' @description
//'
//' Kernel functions that transform observed p-values or their support according
//' to a discrete FWER approach. The outputs are used by [`discrete.FWER()`].
//' `kernel_DFWER_fast`, computes the transformed \eqn{p}-values, while
//' `kernel_DFWER_crit` additionally computes and returns the critical
//' constants. The end user should not use these functions directly, as they are
//' internal functions and parameters (including their names, order, etc.) may
//' be changed without notice!
//' 
//' @templateVar pCDFlist TRUE
//' @template param
//' 
//' @param pvalues        numeric vector, sorted in increasing order, that
//'                       either must contain the entirety of all observable
//'                       values of the p-value supports (when computing
//'                       critical constants) or only the sorted raw p-values.
//' @param independence   single boolean specifying whether the \eqn{p}-values
//'                       are independent; if FALSE (the default), the discrete
//'                       Bonferroni procedure \[d-Bonf\] is performed;
//'                       otherwise, \[d-Ind\] is computed.
//' @param pCDFcounts     integer vector of counts that indicates to how many
//'                       p-values each **unique** p-value distributions
//'                       belongs.
//' @param support        numeric vector, sorted in increasing order, that
//'                       contains the entirety of all observable values of the
//'                       p-value supports.
//' @param sorted_pv      numeric vector, sorted in increasing order, containing
//'                       the raw p-values.
//' @param alpha          single real number strictly between 0 and 1 indicating
//'                       the target FWER level.
//' @param pCDFindices    list of integer vectors containing the indices that
//'                       indicate to which raw \eqn{p}-value in `sorted_pv`
//'                       each item in `pCDFlist` belongs, and must have the
//'                       same length as `pCDFlist`; if `NULL` (the default), it
//'                       is assumed that the first item of `pCDFlist`
//'                       corresponds to the first \eqn{p}-value, the second
//'                       item to the second \eqn{p}-value etc. in which case
//'                       the lengths of `pCDFlist` and `sorted_pv` must be
//'                       equal.
//' 
//' @return
//' For `kernel_DFWER_singlestep_fast()` and `kernel_DFWER_stepwise_fast()` a
//' vector of transformed p-values is returned. `kernel_DFWER_singlestep_crit`
//' and `kernel_DFWER_stepwise_crit` return a list with critical constants
//' (`$crit_consts`) and adjusted p-values (`$pval_transf`).
//' 
//' @seealso
//' [`discrete_FWER()`], [`direct_discrete_FWER()`]
//' 
//' @template example
//' @examples
//' alpha <- 0.05
//' 
//' //' # If not searching for critical constants, we only get the adjusted p-values
//' sorted_pvals <- sort(raw_pvalues)
//' y_dBonf_fast <- DiscreteFWER:::kernel_DFWER_singlestep_fast(pCDFlist, sorted_pvals)
//' y_dSid_fast  <- DiscreteFWER:::kernel_DFWER_singlestep_fast(pCDFlist, sorted_pvals, TRUE)
//' # Adjusted values
//' y_dBonf_fast
//' y_dSid_fast
//' 
//' # Compute transformed support
//' support      <- unique(sort(unlist(pCDFlist)))
//' y_dBonf_crit <- DiscreteFWER:::kernel_DFWER_crit(pCDFlist, support, sorted_pvals, alpha)
//' y_dSid_crit  <- DiscreteFWER:::kernel_DFWER_crit(pCDFlist, support, sorted_pvals, alpha, TRUE)
//' # critical constants
//' y_dBonf_crit$crit_consts
//' y_dSid_crit$crit_consts
//' # Adjusted p-values
//' y_dBonf_crit$pval_transf
//' y_dSid_crit$pval_transf
//' 
//'

///' @export
//' @rdname kernel
// [[Rcpp::export]]
NumericVector kernel_DFWER_singlestep_fast(const List& pCDFlist, const NumericVector& pvalues, const bool independence = false, const Nullable<IntegerVector>& pCDFcounts = R_NilValue);

///' @export
//' @rdname kernel
// [[Rcpp::export]]
List kernel_DFWER_singlestep_crit(const List& pCDFlist, const NumericVector& support, const NumericVector& sorted_pv, const double alpha = 0.05, const bool independence = false, const Nullable<IntegerVector>& pCDFcounts = R_NilValue);

///' @export
//' @rdname kernel
// [[Rcpp::export]]
NumericVector kernel_DFWER_stepwise_fast(const List& pCDFlist, const NumericVector& sorted_pv, const bool independence = false, const Nullable<List>& pCDFindices = R_NilValue);

///' @export
//' @rdname kernel
// [[Rcpp::export]]
List kernel_DFWER_stepwise_crit(const List& pCDFlist, const NumericVector& support, const NumericVector& sorted_pv, const double alpha = 0.05, const bool independence = false, const Nullable<List>& pCDFindices = R_NilValue);
