#' The transmission matrix for a common genetic model
#'
#' A function to calculate the transmission matrix for a single autosomal
#' genetic locus with an arbitrary number of alleles and unphased genotypes,
#' based on Mendel's laws of inheritance.
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
#' There are `ngeno = n_alleles*(n_alleles + 1)/2` possible unphased genotypes,
#' and by choosing an order on these genotypes (see below) we can take the set
#' of possible genotypes to be `1:ngeno`.
#' Then the `(ngeno * gm + gf - ngeno, go)`th element of the outputted matrix is
#' the conditional probability that a person has genotype `go`, given that his
#' or her biological mother and father have genotypes `gm` and `gf`,
#' respectively.
#' When `annotate` is `TRUE`, the function converts this matrix to a data frame,
#' adds column names giving the offspring genotype corresponding to each
#' column, and adds columns `gm` and `gf` describing the parental genotypes
#' corresponding to each row.
#' In this data frame, all genotypes are written in the usual form
#' `1/1, 1/2, ...` for alleles `1:n_alleles`, so these
#' annotations show the genotype order referred to above.
#' Note that if the output of this function is to be used as the `trans`
#' argument of \code{\link{pedigree_loglikelihood}} then the `annotate` option
#' must be set to `FALSE`.
#'
#' @return Either a matrix of genetic transmission probabilities suitable to be
#' used as the `trans` argument of \code{\link{pedigree_loglikelihood}}
#' (if `annotate` is `FALSE`), or a data frame that is an annotated version of
#' this matrix (if `annotate` is `TRUE`).
#'
#' @export
#'
#' @examples
#' # The transition matrix for a biallelic, autosomal locus with unphased genotypes
#' trans_monogenic(2)
#' trans_monogenic(2, annotate = TRUE)
#'
trans_monogenic <- function(n_alleles, annotate = FALSE) {
  # Calculate the transmission matrix for a single, autosomal genetic locus with
  # n_alleles alleles and unphased genotypes
  # n_alleles <- 3

  # List the possible (unordered) genotypes geno.list
  # The possible alleles are 1:n_alleles and the possible genotypes are
  # 1:(n_alleles*(n_alleles+1)/2),
  # with genotype i corresponds to the unphased genotype geno.list[i]
  # eg <- expand.grid(1:n_alleles,1:n_alleles)
  eg <- cbind(rep(1:n_alleles, each = n_alleles),
              rep(1:n_alleles, times = n_alleles))
  f <- function(x) {
    if (x[1] > x[2]) {
      return("")
    } else {
      return(paste(x, collapse = "/"))
    }
  }
  geno.list <- apply(eg, 1, f)
  geno.list <- geno.list[geno.list != ""]
  ng <- length(geno.list)
  geno.list
  ng - n_alleles * (n_alleles + 1) / 2

  # Calculate the transmission matrix
  trans <- matrix(0, ng^2, ng)
  gm.list <- rep("", ng^2)
  gf.list <- rep("", ng^2)
  for (gm in 1:ng) {
    for (gf in 1:ng) {
      # gm <- 1; gf <- 3
      am.list <- strsplit(geno.list[gm], "/", fixed = TRUE)[[1]]
      af.list <- strsplit(geno.list[gf], "/", fixed = TRUE)[[1]]
      geno.off <- function(am, af) {
        if (am <= af) {
          go <- paste(c(am, af), collapse = "/")
        } else {
          go <- paste(c(af, am), collapse = "/")
        }
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
