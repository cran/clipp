#' clipp: Calculate Likelihoods by Pedigree Paring
#'
#' `clipp` provides a fast and general implementation of the Elston-Stewart
#' algorithm, and can calculate the log-likelihoods of large and complex pedigrees
#' without loops.  General references for the Elston-Stewart algorithm are
#' (Elston & Stewart, 1971), (Lange & Elston, 1975) and (Cannings et al., 1978).
#'
#' The main function is \code{\link{pedigree_loglikelihood}}, which calculates
#' the pedigree likelihood on page 117 of (Lange, 2002) for almost any choice of
#' genotype frequencies, transmission matrix and penetrance matrix.  Helper
#' functions are provided to calculate the genotype frequencies and transmission
#' matrices for genetic models that often arise in applications.  The function
#' \code{\link{genotype_probabilities}} calculates genotype probabilities
#' for a target person within a family, given the family's phenotypes.
#'
#' The current implementation of `clipp` does not allow pedigree loops, such as
#' those caused by inbreeding or by two sisters having children with two brothers
#' from a different family (see (Totir et al., 2009) for a precise definition).
#'
#' It is feasible to apply `clipp` to very large families,
#' e.g. in the examples for \code{\link{pedigree_loglikelihood}},
#' the log-likelihood of one family with approximately 10,000 members is calculated
#' in less than one minute on a standard desktop computer.
#' Numerical issues will eventually limit the family size,
#' though `clipp` takes care to avoid arithmetic underflow and other issues.
#'
#' @references
#' Cannings C, Thompson E, Skolnick M. Probability functions
#' on complex pedigrees. Advances in Applied Probability, 1978;10(1):26-61.
#'
#' Elston RC, Stewart J. A general model for the genetic analysis of pedigree
#' data. Hum Hered. 1971;21(6):523-542.
#'
#' Lange K.  Mathematical and Statistical Methods for Genetic Analysis
#' (second edition). Springer, New York. 2002.
#'
#' Lange K, Elston RC. Extensions to pedigree analysis I. Likehood calculations
#' for simple and complex pedigrees. Hum Hered. 1975;25(2):95-105.
#'
#' Totir LR, Fernando RL, Abraham J. An efficient algorithm to compute marginal
#' posterior genotype probabilities for every member of a pedigree with loops.
#' Genet Sel Evol. 2009;41(1):52.
#'
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
