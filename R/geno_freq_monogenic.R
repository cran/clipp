#' Calculate genotype frequencies from allele frequencies using Hardy-Weinberg equilibrium
#'
#' A function to calculate the unphased genotype frequencies for a single
#' autosomal genetic locus that has given allele frequencies and is at
#' Hardy-Weinberg equilibrium (HWE).
#'
#' @param p_alleles  A vector of strictly positive numbers that sum to `1`, with
#' `p_alleles[i]` interpreted as the allele frequency of the `i`th allele of the
#' genetic locus.  When `annotate` is `TRUE`, the names of the alleles will be
#' taken to be `names(p_alleles)` or, if `names(p_alleles)` is `NULL`, to be
#' `1:length(p_alleles)`.
#'
#' @param annotate  A logical flag.  When `FALSE` (the default), the function
#' returns a vector suitable to be used as the `geno_freq` argument of
#' \code{\link{pedigree_loglikelihood}}.  When `TRUE`, the function adds a
#' `names` attribute to this vector to indicate which genotype corresponds
#' to which element.
#'
#' @details  For a genetic locus at Hardy-Weinberg equilibrium, the population
#' allele frequencies at the locus determine the population genotype frequencies;
#' see Sections 1.2 and 1.3 of (Lange, 2002).
#' Given a vector `p_alleles` containing the allele
#' frequencies, this function returns the frequencies of the possible unphased
#' genotypes, in a particular order that can be viewed by setting `annotate` to
#' `TRUE`. If the alleles are named `1:length(p_alleles)`, so that `p_alleles[i]`
#' is the frequency of allele `i`, then the unphased genotypes are named
#' `1/1, 1/2, ...`.  Note that if the output of this function is to be used
#' as the `geno_freq` argument of \code{\link{pedigree_loglikelihood}}
#' then the `annotate` option must be set to `FALSE`.
#'
#' @return A vector of strictly positive numbers (the genotype frequencies)
#' that sum to `1`, named with the genotype names if `annotate` is `TRUE`.
#'
#' @export
#'
#' @examples
#' # Genotype frequencies for a biallelic locus at HWE and with a minor allele frequency of 10%
#' p_alleles <- c(0.9, 0.1)
#' geno_freq_monogenic(p_alleles, annotate = TRUE)
#'
#' # Genotype frequencies for a triallelic locus at Hardy-Weinberg equilibrium
#' p_alleles <- c(0.85, 0.1, 0.05)
#' geno_freq_monogenic(p_alleles, annotate = TRUE)
#' sum(geno_freq_monogenic(p_alleles))
#'
#' @references
#' Lange K.  Mathematical and Statistical Methods for Genetic Analysis (second edition).
#' Springer, New York.  2002.
#'
geno_freq_monogenic <- function(p_alleles, annotate = FALSE) {
  # The genotype frequencies corresponding to the allele frequencies p_alleles under Hardy-Weinberg equilibrium
  # The elements of p_alleles must sum to 1 and be strictly positive

  # List the possible (unordered) genotypes geno.list, as in trans.unphased()
  n_alleles <- length(p_alleles)
  # eg <- expand.grid(1:n_alleles,1:n_alleles)
  eg <- cbind(rep(1:n_alleles, each = n_alleles), rep(1:n_alleles, times = n_alleles))
  f1 <- function(x) {
    if (x[1] > x[2]) {
      return("")
    } else {
      return(paste(x, collapse = "/"))
    }
  }
  geno.list <- apply(eg, 1, f1)
  geno.list <- geno.list[geno.list != ""]
  ng <- length(geno.list)

  # Calculate the genotype probabilities
  pg <- numeric(ng)
  for (j in 1:ng) {
    a.list <- as.numeric(strsplit(geno.list[j], "/", fixed = T)[[1]])
    if (a.list[1] == a.list[2]) {
      pg[j] <- p_alleles[a.list[1]] * p_alleles[a.list[2]]
    } else {
      pg[j] <- 2 * p_alleles[a.list[1]] * p_alleles[a.list[2]]
    }
  }

  # Add genotype names to the list?
  allele.names <- 1:n_alleles
  if (!is.null(names(p_alleles)))  {allele.names <- names(p_alleles)}
  f2 <- function(j) {
    a.list <- as.numeric(strsplit(geno.list[j], "/", fixed = T)[[1]])
    paste0(allele.names[a.list[1]], "/", allele.names[a.list[2]])
  }
  geno.list.names <- sapply(1:ng, f2)
  if (annotate) names(pg) <- geno.list.names

  return(pg)
}
