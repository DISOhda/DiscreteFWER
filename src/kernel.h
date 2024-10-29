#include "helper.h"

//' @name kernel
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
//' @param pvalues       numeric vector, sorted in increasing order, that either
//'                      must contain the entirety of all observable values of
//'                      the p-value supports (when computing critical
//'                      constants) or only the sorted raw p-values.
//' @param independent   single boolean specifying whether the \eqn{p}-values
//'                      are independent; if FALSE (the default), the discrete
//'                      Bonferroni procedure \[d-Bonf\] is performed;
//'                      otherwise, \[d-Ind\] is computed.
//' @param pCDFcounts    integer vector of counts that indicates to how many
//'                      p-values each **unique** p-value distributions belongs.
//' @param support       numeric vector, sorted in increasing order, that
//'                      contains the entirety of all observable values of the
//'                      p-value supports.
//' @param sorted_pv     numeric vector, sorted in increasing order, containing
//'                      the raw p-values.
//' @param alpha         single real number strictly between 0 and 1 indicating
//'                      the target FWER level.
//' 
//' @return
//' For `kernel_DFWER_fast()` a vector of transformed p-values is returned.
//' `kernel_DFWER_crit` returns a list with critical constants (`$crit.consts`)
//' and transformed p-values (`$pval.transf`).
//' 
//' @seealso
//' [`discrete.FWER()`], [`direct.discrete.FWER()`]
//' 
//' @template example
//' @examples
//' alpha <- 0.05
//' 
//' \dontrun{
//' # If not searching for critical constants, we use only the observed p-values
//' sorted.pvals <- sort(raw.pvalues)
//' y.dBonf.fast <- DiscreteFWER:::kernel_DFWER_fast(pCDFlist, sorted.pvals)
//' y.dInd.fast  <- DiscreteFWER:::kernel_DFWER_fast(pCDFlist, sorted.pvals, TRUE)
//' # transformed values
//' y.dBonf.fast
//' y.dInd.fast
//' 
//' # compute transformed support
//' support      <- sort(unique(unlist(pCDFlist)))
//' y.dBonf.crit <- DiscreteFWER:::kernel_DFWER_crit(pCDFlist, support, sorted.pvals, alpha)
//' y.dInd.crit  <- DiscreteFWER:::kernel_DFWER_crit(pCDFlist, support, sorted.pvals, alpha, TRUE)
//' # critical constants
//' y.dBonf.crit$crit.consts
//' y.dInd.crit$crit.consts
//' # Transformed p-values
//' y.dBonf.crit$pval.transf
//' y.dInd.crit$pval.transf
//' }
//'

///' @export
//' @rdname kernel
// [[Rcpp::export]]
NumericVector kernel_DFWER_single_fast(const List &pCDFlist, const NumericVector &pvalues, const bool independent = false, const Nullable<NumericVector> &pCDFcounts = R_NilValue);

///' @export
//' @rdname kernel
// [[Rcpp::export]]
List kernel_DFWER_single_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const double alpha = 0.05, const bool independent = false, const Nullable<NumericVector> &pCDFcounts = R_NilValue);
