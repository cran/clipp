#' Genotype frequencies for the hypergeometric polygenic model
#'
#' A function to calculate the genotype frequencies for the hypergeometric
#' polygenic model of (Cannings et al., 1978), see Section 8.9 of (Lange, 2002)
#' for a nice description of this model.
#'
#' @param n_loci A positive integer, interpreted as the number of biallelic
#' genetic loci that contribute to the polygene.  The polygene will have
#' `2*n_loci + 1` genotypes, so `n_loci` is typically fairly small, e.g. `4`.
#'
#' @param annotate  A logical flag.  When `FALSE` (the default), the function
#' returns a vector suitable to be used as the `geno_freq` argument of
#' \code{\link{pedigree_loglikelihood}}.  When `TRUE`, the function adds a
#' `names` attribute to this vector to indicate which element corresponds to
#' which genotype.
#'
#' @details  The hypergeometric polygenic model
#' (Cannings et al., 1978; Lange, 2002) is a computationally feasible
#' genetic model that approximates the combined effect of a given number
#' (`n_loci`) of unlinked biallelic genetic loci.  This model is often used to model
#' the effect of such loci on a trait when the alleles at these loci either
#' increase the trait by a certain, locus-independent amount (if a 'positive'
#' allele) or decrease
#' the trait by the same amount (if a 'negative' allele), with 'positive' and
#' 'negative' alleles equally likely at each locus.  In this case, the only
#' relevant aspect of the `3 ^ n_loci` possible joint genotypes is the total number
#' of 'positive' alleles, so the possible genotypes of
#' the hypergeometric polygenic model are taken to be `0:(2*n_loci)`.
#' The transmission probabilities and genotype frequencies of the hypergeometric
#' polygenic model approximate these quantities for the combination of the `n_loci`
#' biallelic genetic loci described above.  Under this model, the polygenic
#' genotype for each person is approximately normally distributed, and these
#' genotypes are correlated within families with correlation coefficients
#' (in non-inbred families) equal to the kinship coefficients (Lange, 2002).
#'
#' Setting `annotate` to `TRUE` names each element of the output vector with
#' the corresponding genotype. The `annotate` option must be set to `FALSE`
#' if the output of this function is to be used as the `geno_freq` argument of
#' \code{\link{pedigree_loglikelihood}}.
#'
#' @return A vector of strictly positive numbers (the genotype frequencies)
#' that sum to `1`, named with the genotype names if `annotate` is `TRUE`.
#'
#' @examples
#' geno_freq_polygenic(4, annotate = TRUE)
#' sum(geno_freq_polygenic(4))
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
geno_freq_polygenic <- function(n_loci, annotate = FALSE) {
  #n_loci <- 4; annotate <- T

  # List the possible genotypes
  geno.list <- 0:(2*n_loci)

  # Calculate the genotype probabilities
  f <- function(i)  choose(2*n_loci, i) * 2^(-2*n_loci)
  pg <- sapply(geno.list, f)

  # Add genotype names to the list?
  if (annotate)  names(pg) <- geno.list

  return(pg)
}
