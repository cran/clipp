#' The transmission matrix for phased genotypes
#'
#' A function to calculate the transmission matrix for a single autosomal
#' genetic locus with an arbitrary number of alleles and phased genotypes,
#' based on Mendel's laws of inheritance.  Phased genotypes can be used to
#' investigate parent-of-origin effects, e.g. see (van Vliet et al., 2011).
#'
#' @param n_alleles A positive integer, interpreted as the number of possible
#' alleles at the genetic locus.
#'
#' @param annotate A logical flag. When `FALSE` (the default), the function
#' returns a matrix suitable to be used as the `trans` argument of
#' \code{\link{pedigree_loglikelihood}}. When `TRUE`, the function annotates
#' this matrix (and converts it to a data frame) to make the output more
#' easily understood by humans.
#'
#' @details When `annotate` is `FALSE`, a matrix of genetic transmission
#' probabilities is returned, with rows corresponding to the possible joint
#' parental genotypes and columns corresponding to the possible offspring
#' genotypes.
#' There are `ngeno = n_alleles^2` possible phased genotypes,
#' and by choosing an order on these genotypes (which can be
#' viewed by setting `annotate` to `TRUE`, see below)
#' we can label the set of possible phased genotypes as `1:ngeno`.
#' Then the `(ngeno * gm + gf - ngeno, go)`th element of the outputted matrix is
#' the conditional probability that a person has genotype `go`, given that his
#' or her biological mother and father have genotypes `gm` and `gf`,
#' respectively.
#'
#' When `annotate` is `TRUE`, the function converts this matrix to a data frame,
#' adds column names giving the offspring genotype corresponding to each
#' column, and adds columns `gm` and `gf` describing the parental genotypes
#' corresponding to each row.
#' In this data frame, phased genotypes are written in the usual form
#' `1|1, 1|2, ...` for the alleles `1:n_alleles`, where `a|b` means the
#' maternal allele is `a` and the paternal allele is `b`.
#'
#' Note that if the output of this function is to be used as the `trans`
#' argument of \code{\link{pedigree_loglikelihood}} then the `annotate` option
#' must be set to `FALSE`.
#'
#' @return Either a matrix of genetic transmission probabilities suitable to be
#' used as the `trans` argument of \code{\link{pedigree_loglikelihood}}
#' (if `annotate` is `FALSE`), or a data frame that is an annotated version of
#' this matrix (if `annotate` is `TRUE`).
#'
#' @examples
#' # The transition matrix for a biallelic, autosomal locus with phased genotypes
#' trans_phased(2)
#' trans_phased(2, annotate = TRUE)
#'
#' @references
#' van Vliet CM, Dowty JG, van Vliet JL, et al. Dependence of colorectal cancer
#' risk on the parent-of-origin of mutations in DNA mismatch repair genes.
#' Hum Mutat. 2011;32(2):207-212.
#'
#' @export
#'
trans_phased <- function(n_alleles, annotate = FALSE) {
  # Calculate the transmission matrix for a single, autosomal genetic locus with
  # n_alleles alleles and phased genotypes
  # n_alleles <- 2; annotate <- T

  # List the possible phased genotypes geno.list
  # The possible alleles are 1:n_alleles and the possible phased genotypes are
  # 1:(n_alleles^2), with genotype i corresponds to the phased genotype geno.list[i]
  # eg <- expand.grid(1:n_alleles,1:n_alleles)
  eg <- cbind(rep(1:n_alleles, each = n_alleles),
              rep(1:n_alleles, times = n_alleles))
  f <- function(x) {  return(paste(x, collapse = "|"))  }
  geno.list <- apply(eg, 1, f)
  ng <- length(geno.list)
  geno.list
  ng - n_alleles^2

  # Calculate the transmission matrix
  trans <- matrix(0, ng^2, ng)
  gm.list <- rep("", ng^2)
  gf.list <- rep("", ng^2)
  for (gm in 1:ng) {
    for (gf in 1:ng) {
      # gm <- 1; gf <- 3
      am.list <- strsplit(geno.list[gm], "|", fixed = TRUE)[[1]]
      af.list <- strsplit(geno.list[gf], "|", fixed = TRUE)[[1]]
      geno.off <- function(am, af) {
        go <- paste(c(am, af), collapse = "|")
        return(which(geno.list == go))
      }
      for (i in 1:2) {
        for (j in 1:2) {
          go <- geno.off(am.list[i], af.list[j])
          trans[ng * gm + gf - ng, go] <- trans[ng * gm + gf - ng, go] + 0.25
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
