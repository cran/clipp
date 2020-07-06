#' Calculate the genotype probabilities of a target person
#'
#' For a given individual within a given family, calculate the person's
#' conditional genotype probabilities, given the family's phenotypes and
#' relationship structure
#'
#' @param target_ID  The individual identifier (`indiv`) of the person in
#' the pedigree `fam` whose genotype probabilities are sought.
#'
#' @param fam A data frame with rows corresponding to people and columns
#' corresponding to the following variables (other variables can be included
#' but will be ignored), which will be coerced to `character` type:
#'
#' * `indiv`, an identifier for each individual person, with no duplicates
#'    in `fam`.
#' * `mother`, the individual identifier of each person's mother, or missing
#'    (`NA`) for founders.
#' * `father`, the individual identifier of each person's father, or missing
#'    (`NA`) for founders.
#'
#' @param geno_freq  A vector of strictly positive numbers that sum to `1`.
#' The possible genotypes of the underlying genetic model are
#' `1:length(geno_freq)`, and `geno_freq[j]` is interpreted as the population
#' frequency of genotype `j`. For certain genetic models that often
#' occur in applications, these genotype frequencies can be calculated by
#' \code{\link{geno_freq_HWE}}.
#'
#' @param trans  An `ngeno^2` by `ngeno` matrix of non-negative numbers whose rows
#' all sum to `1`, where `ngeno = length(geno_freq)` is the number of possible
#' genotypes.  The rows of `trans` correspond to joint parental
#' genotypes and the columns to offspring genotypes.
#' The element `trans[ngeno * gm + gf - ngeno, go]` is interpreted as the conditional
#' probability that a person has genotype `go`, given that his or her
#' biological mother and father have genotypes `gm` and `gf`, respectively.
#' For certain genetic models that often occur in applications,
#' this transmission matrix can be calculated by \code{\link{trans_monogenic}}.
#'
#' @param penet An `nrow(dat)` by `length(geno_freq)` matrix of non-negative
#' numbers. The element `penet[i,j]` is interpreted as the conditional
#' probability of the phenotypes of the person corresponding to row `i` of
#' `fam`, given that his or her genotype is `j`.
#' Note that known genotype data can be incorporated into `penet` by regarding
#' known genotypes as phenotypes, i.e. by regarding known genotypes as
#' (possibly noisy) measurements of the underlying true genotypes.
#' If any row of `penet` consists entirely of zeroes then the family's
#' phenotypes are impossible, and `-Inf` is returned.
#'
#' @details
#' The genotype probabilities are calculated by essentially the same algorithm
#' as the one that underlies \code{\link{pedigree_loglikelihood}};
#' see there for details.
#'
#' @return A vector of length `length(geno_freq)` whose `j`th element is
#' the conditional probability that the target person has genotype `j`,
#' given the family's relationship structure and phenotypes.
#'
#' @export
#'
#' @examples
#' # Read in some sample data
#' data("dat_small", "penet_small")
#' str(dat_small)
#' str(penet_small)
#'
#' # Calculate the genotype probabilities for individual "ora008" in the family "ora"
#' w <- which(dat_small$family == "ora")
#' fam <- dat_small[w, -1]
#' penet <- penet_small[w, ]
#' trans <- trans_monogenic(2)
#' geno_freq <- geno_freq_HWE(p_alleles = c(0.9, 0.1))
#' genotype_probabilities(target_ID = "ora008", fam, geno_freq, trans, penet)
#'
genotype_probabilities <- function(target_ID, fam, geno_freq, trans, penet) {

  dat <- cbind(family = 1, fam)

  if (!target_ID %in% dat$indiv) {
    stop("Target ID cannot be found in the family data",  call. = FALSE)
  }

  # This saved which row corresponds to target_ID
  # Because target_ID will be converted to numeric later
  target_id <- which(dat$indiv == target_ID)

  # Convert IDs
  dat$mother[which(is.na(dat$mother))] <- ""
  dat$father[which(is.na(dat$father))] <- ""
  dat <- convert_IDs(dat, convert.IDs.numeric = TRUE)

  fam_penet <- list(fam = dat, penet = penet)

  pedigree_loglikelihood_g(
    fam_penet,
    geno_freq,
    trans = trans,
    cl = NULL,
    target_id
  )


}

