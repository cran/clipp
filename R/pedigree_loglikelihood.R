#' Calculate the log-likelihoods of pedigrees
#'
#' For one or more pedigrees, this function calculates the natural logarithm of
#' the pedigree likelihood on page 117 of (Lange, 2002), given inputs that
#' correspond to the terms in this formula.
#'
#' @param dat A data frame with rows corresponding to people and columns
#' corresponding to the following variables (other variables can be included
#' but will be ignored), which will be coerced to `character` type:
#'
#' * `family` (optional), an identifier for each person's family, constant
#'    within families.  If this variable is not supplied then `dat` will be
#'    treated as a single pedigree.
#' * `indiv`, an identifier for each individual person, with no duplicates
#'    across the dataset.
#' * `mother`, the individual identifier of each person's mother, or missing
#'    (`NA`) for founders.
#' * `father`, the individual identifier of each person's father, or missing
#'    (`NA`) for founders.
#'
#' @param geno_freq  A vector of strictly positive numbers that sum to `1`.
#' The possible genotypes of the underlying genetic model are
#' `1:length(geno_freq)`, and `geno_freq[j]` is interpreted as the population
#' frequency of genotype `j`, so `geno_freq` is essentially the
#' function `Prior` of (Lange, 2002). For certain genetic models that often
#' occur in applications, these genotype frequencies can be calculated by
#' \code{\link{geno_freq_HWE}}.
#'
#' @param trans  An `ngeno^2` by `ngeno` matrix of non-negative numbers whose rows
#' all sum to `1`, where `ngeno = length(geno_freq)` is the number of possible
#' genotypes. The rows of `trans` correspond to joint parental genotypes and
#' the columns to offspring genotypes.  The element
#' `trans[ngeno * gm + gf - ngeno, go]` is interpreted as the conditional
#' probability that a person has genotype `go`, given that his or her
#' biological mother and father have genotypes `gm` and `gf`, respectively.
#' So `trans` is essentially the transmission function `Tran` of (Lange, 2002).
#' For certain genetic models that often occur in applications,
#' this transmission matrix can be calculated by \code{\link{trans_monogenic}}.
#'
#' @param penet An `nrow(dat)` by `length(geno_freq)` matrix of non-negative
#' numbers. The element `penet[i,j]` is interpreted as the conditional
#' probability of the phenotypes of the person corresponding to row `i` of
#' `dat`, given that his or her genotype is `j`.  So `penet` is essentially the
#' penetrance function `Pen` of (Lange, 2002).
#' Note that genotype data can be incorporated into `penet` by regarding
#' known genotypes as phenotypes, i.e. by regarding known genotypes as
#' (possibly noisy) measurements of the underlying true genotypes.
#' If any row of `penet` consists entirely of zeroes then the likelihood is `0`,
#' so the returned log-likelihood is `-Inf`.
#'
#' @param sum_loglik A logical flag.  Return a named vector giving the
#' log-likelihood of each family if `sum_loglik` is `FALSE`, or return the sum
#' of these log-likelihoods if `sum_loglik` is `TRUE` (the default).
#'
#' @param ncores The number of cores to be used, with `ncores = 1` (the
#' default) corresponding to non-parallel computing.  When `ncores > 1`,
#' the `parallel` package is used to parallelize the calculation by dividing
#' the pedigrees among the different cores.
#'
#' @details
#' This function provides a fast and general implementation of the
#' Elston-Stewart algorithm to calculate the log-likelihoods of large and
#' complex pedigrees without loops.  General references for the Elston-Stewart
#' algorithm are (Elston & Stewart, 1971), (Lange & Elston, 1975) and
#' (Cannings et al., 1978).
#'
#' Each family within `dat` should be a complete pedigree, meaning that each
#' non-missing mother or father ID should correspond to a row, and each person
#' should either have both parental IDs missing (if a founder) or non-missing
#' (if a non-founder).  No family should contain pedigree loops, such as those
#' caused by inbreeding or by two sisters having children with two brothers
#' from a different family; see (Totir et al., 2009) for a precise definition.
#' The code is not currently suitable for families containing identical twins,
#' since these are not distinguished from full siblings.
#'
#' In `geno_freq`, `trans` and `penet`, the order of the possible genotypes
#' must match (so that the same genotype corresponds to element `j` of `geno_freq`
#' and column `j` of `trans` and `penet`, for each `j` in `1:length(geno_freq)`).
#'
#' @return Either a named vector giving the log-likelihood of each family
#' or the sum of these log-likelihoods, depending on `sum_loglik` (see above).
#'
#' @export
#'
#' @examples
#' # Load pedigree files and penetrance matrices
#' data("dat_small", "penet_small", "dat_large", "penet_large")
#'
#' # Settings for a biallelic locus in Hardy-Weinberg equilibrium and with a
#' # minor allele frequency of 10%
#' geno_freq <- geno_freq_HWE(c(0.9, 0.1))
#' trans <- trans_monogenic(2)
#'
#' # Calculate the log-likelihoods for 10 families, each with approximately 100 family members
#' pedigree_loglikelihood(
#'   dat_small, geno_freq, trans, penet_small, sum_loglik = FALSE, ncores = 2
#' )
#'
#' # Calculate the log-likelihood for one family with approximately 10,000 family members
#' # Note:  this calculation takes approximately one minute on a standard desktop computer
#' # Note:  parallelization achieves nothing here because there is only one family
#' str(dat_large)
#' \donttest{
#' system.time(
#'   ll <- pedigree_loglikelihood(dat_large, geno_freq, trans, penet_large)
#' )
#' ll
#' }
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
pedigree_loglikelihood <- function(dat, geno_freq, trans, penet,
                                   sum_loglik = TRUE, ncores = 1) {

  #check_dat(dat)

  # Check if family can be found
  one_fam <- FALSE
  if (!"family" %in% names(dat)) {
    one_fam <- TRUE
    dat <- cbind(family = 1, dat)
  }

  # Convert IDs
  dat$mother[which(is.na(dat$mother))] <- ""
  dat$father[which(is.na(dat$father))] <- ""
  dat <- convert_IDs(dat, convert.IDs.numeric = TRUE)

  # Prepare a list of fam and penet for parallelization
  ufamID <- unique(dat$family)
  nfam <- length(ufamID)
  fam_penet_list <- vector("list", length = nfam)
  for (i in seq_along(ufamID)) {
    famID <- ufamID[i]
    keep <- which(dat$family == famID)
    fam <- dat[keep, -1]
    peneti <- penet[keep, ]
    fam_penet_list[[i]] <- list(
      fam = fam,
      penet = peneti
    )
  }
  names(fam_penet_list) <- ufamID

  if (ncores == 1) {
    ll <- lapply(fam_penet_list, pedigree_loglikelihood_g,
                 geno_freq = geno_freq, trans = trans)
  } else {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    ll <- parallel::parLapply(
      cl = cl, X = fam_penet_list, fun = pedigree_loglikelihood_g,
      geno_freq = geno_freq, trans = trans)
  }

  ll <- unlist(ll)
  if (sum_loglik) {
    sum(ll)
  } else {
    #df <- data.frame(family = ufamID, loglikelihood = ll)
    #rownames(df) <- NULL
    if(one_fam) {
      names(ll) <- NULL
    } else {
      names(ll) <- ufamID
    }
    ll
  }
}
