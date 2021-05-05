#' The transmission matrix for the hypergeometric polygenic model
#'
#' A function to calculate the transmission matrix for the hypergeometric
#' polygenic model of (Cannings et al., 1978), see also Section 8.9 of
#' (Lange, 2002) for a nice description of this model.
#'
#' @param n_loci A positive integer, interpreted as the number of biallelic
#' genetic loci that contribute to the polygene.  The polygene will have
#' `2*n_loci + 1` genotypes, so `n_loci` is typically fairly small, e.g. `4`.
#'
#' @param annotate A logical flag. When `FALSE` (the default), the function
#' returns a matrix suitable to be used as the `trans` argument of
#' \code{\link{pedigree_loglikelihood}}. When `TRUE`, the function annotates
#' this matrix (and converts it to a data frame) to make the output more
#' easily understood by humans.
#'
#' @details This function calculates the genetic transmission probabilities
#' (i.e. the conditional probability of a person's genotype, given his
#' or her biological parents' genotypes) for the hypergeometric polygenic model,
#' which is described in \code{\link{geno_freq_polygenic}}.
#'
#' When `annotate` is `FALSE`,
#' a matrix of transmission probabilities is returned, with rows
#' corresponding to the possible joint parental genotypes and columns
#' corresponding to the possible offspring genotypes.
#' Setting `annotate` to `TRUE` shows which rows and columns correspond to
#' which genotypes,  by adding offspring genotypes as column names and adding
#' columns `gm` and `gf` containing (respectively) the mother's and father's
#' genotypes. Note that if the output of this function is to be used as the `trans`
#' argument of \code{\link{pedigree_loglikelihood}} then the `annotate` option
#' must be set to `FALSE`.
#'
#' @return Either a matrix of genetic transmission probabilities suitable to be
#' used as the `trans` argument of \code{\link{pedigree_loglikelihood}}
#' (if `annotate` is `FALSE`), or a data frame that is an annotated version of
#' this matrix (if `annotate` is `TRUE`).
#'
#' @examples
#' trans_polygenic(4, annotate = TRUE)
#' apply(trans_polygenic(4), 1, sum)
#'
#' @references
#' Cannings C, Thompson E, Skolnick M. Probability functions
#' on complex pedigrees. Advances in Applied Probability, 1978;10(1):26-61.
#'
#' Lange K.  Mathematical and Statistical Methods for Genetic Analysis
#' (second edition). Springer, New York. 2002.
#'
#' @export
#'
trans_polygenic <- function(n_loci, annotate = FALSE) {
  # n_loci <- 4; annotate <- T

  # List the possible genotypes
  ng <- 2*n_loci + 1
  geno.list <- 0:(ng-1)

  # Calculate the transmission matrix
  f <- function(i,j)  {choose(i,j) * choose(2*n_loci-i, n_loci-j) / choose(2*n_loci, n_loci)}
  trans <- matrix(0, ng^2, ng)
  gm.list <- rep("", ng^2)
  gf.list <- rep("", ng^2)
  for (gm in 1:ng) {
    for (gf in 1:ng) {
      for (hm in 0:n_loci) {
        for (hf in 0:n_loci) {
          p <- f(gm-1, hm) * f(gf-1, hf)
          go <- hm + hf + 1
          trans[ng * gm + gf - ng, go] <- trans[ng * gm + gf - ng, go] + p
        }
      }
      gm.list[ng * gm + gf - ng] <- geno.list[gm]
      gf.list[ng * gm + gf - ng] <- geno.list[gf]
    }
  }

  if (annotate) {
    trans <- data.frame(trans)
    trans <- data.frame(gm = gm.list, gf = gf.list, trans)
    names(trans)[3:ncol(trans)] <- geno.list
  }

  return(trans)
}
