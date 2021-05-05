#' Calculate phased genotype frequencies from allele frequencies,
#' assuming Hardy-Weinberg equilibrium
#'
#' A function to calculate the population frequencies of the phased genotypes
#' at a single autosomal genetic locus that has given allele frequencies and
#' is at Hardy-Weinberg equilibrium.  Phased genotypes can be used to
#' investigate parent-of-origin effects, e.g. see (van Vliet et al., 2011).
#'
#' @param p_alleles  A vector of strictly positive numbers that sum to `1`,
#' with `p_alleles[i]` interpreted as the allele frequency of the `i`th allele
#' of the genetic locus.  When `annotate` is `TRUE`, the names of the alleles
#' will be taken to be `names(p_alleles)` or, if `names(p_alleles)` is `NULL`,
#' to be `1:length(p_alleles)`.
#'
#' @param annotate  A logical flag.  When `FALSE` (the default), the function
#' returns a vector suitable to be used as the `geno_freq` argument of
#' \code{\link{pedigree_loglikelihood}}.  When `TRUE`, the function adds a
#' `names` attribute to this vector to indicate which element corresponds to
#' which phased genotype.
#'
#' @details  For a genetic locus that is at Hardy-Weinberg equilibrium in a
#' particular population, the population allele frequencies at
#' the locus determine the population genotype frequencies; see Sections 1.2 and
#' 1.3 of (Lange, 2002) for the unphased version of this law.  When a genetic
#' locus is at Hardy-Weinberg equilibrium, the maternal and paternal alleles of a
#' random person from the population are independent.  A phased genotype at a
#' genetic locus is an ordered pair consisting of a maternal
#' and paternal allele at the locus.  So to any heterozygous unphased genotype,
#' there are two corresponding phased genotypes, and these two phased genotypes
#' have equal frequencies under Hardy-Weinberg equilibrium.
#'
#' Given a vector `p_alleles` containing the allele frequencies,
#' this function returns the frequencies of the possible phased genotypes,
#' in a particular order that can be viewed by setting `annotate` to `TRUE`.
#' If the alleles are named `1:length(p_alleles)`, so that `p_alleles[i]` is
#' the frequency of allele `i`, then the phased genotypes are of the form
#' `1|1, 1|2, ...`, where `a|b` means the maternal allele is `a` and the
#' paternal allele is `b`. Note that if the output of this function is to be
#' used as the `geno_freq` argument of \code{\link{pedigree_loglikelihood}}
#' then the `annotate` option must be set to `FALSE`.
#'
#' @return A vector of strictly positive numbers (the genotype frequencies)
#' that sum to `1`, named with the genotype names if `annotate` is `TRUE`.
#'
#' @examples
#' # Genotype frequencies for a biallelic locus at Hardy-Weinberg equilibrium
#' # and with a minor allele frequency of 10%
#' p_alleles <- c(0.9, 0.1)
#' geno_freq_phased(p_alleles, annotate = TRUE)
#'
#' # Genotype frequencies for a triallelic locus at Hardy-Weinberg equilibrium
#' p_alleles <- c(0.85, 0.1, 0.05)
#' geno_freq_phased(p_alleles, annotate = TRUE)
#' sum(geno_freq_phased(p_alleles))
#'
#' @references
#' Lange K.  Mathematical and Statistical Methods for Genetic Analysis (second edition).
#' Springer, New York.  2002.
#'
#' van Vliet CM, Dowty JG, van Vliet JL, et al. Dependence of colorectal cancer
#' risk on the parent-of-origin of mutations in DNA mismatch repair genes.
#' Hum Mutat. 2011;32(2):207-212.
#'
#' @export
#'
geno_freq_phased <- function(p_alleles, annotate = FALSE) {
  # The genotype frequencies corresponding to the allele frequencies p_alleles
  # under Hardy-Weinberg equilibrium.  The elements of p_alleles must sum to 1
  # and be strictly positive
  #p_alleles <- c(0.9, 0.07, 0.03); annotate <- T

  # List the possible (unordered) genotypes geno.list, as in trans.unphased()
  n_alleles <- length(p_alleles)
  # eg <- expand.grid(1:n_alleles,1:n_alleles)
  eg <- cbind(rep(1:n_alleles, each = n_alleles), rep(1:n_alleles, times = n_alleles))
  f1 <- function(x) {  return(paste(x, collapse = "|"))  }
  geno.list <- apply(eg, 1, f1)
  ng <- length(geno.list)

  # Calculate the genotype probabilities
  pg <- numeric(ng)
  for (j in 1:ng) {
    a.list <- as.numeric(strsplit(geno.list[j], "|", fixed = T)[[1]])
    pg[j] <- p_alleles[a.list[1]] * p_alleles[a.list[2]]
  }

  # Add genotype names to the list?
  allele.names <- 1:n_alleles
  if (!is.null(names(p_alleles)))  {allele.names <- names(p_alleles)}
  f2 <- function(j) {
    a.list <- as.numeric(strsplit(geno.list[j], "|", fixed = T)[[1]])
    paste0(allele.names[a.list[1]], "|", allele.names[a.list[2]])
  }
  geno.list.names <- sapply(1:ng, f2)
  if (annotate) names(pg) <- geno.list.names

  return(pg)
}
