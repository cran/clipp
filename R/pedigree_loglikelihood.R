#' Calculate the log-likelihoods of pedigrees
#'
#' For one or more pedigrees, this function calculates the natural logarithm of
#' the pedigree likelihood that is on page 117 of (Lange, 2002), given inputs
#' that correspond to the terms in this formula.
#'
#' @param dat A data frame with rows corresponding to people and columns
#' corresponding to the following variables (other variables can be included
#' but will be ignored), which will be coerced to `character` type:
#'
#' * `family` (optional), an identifier for each person's family, constant
#'    within families.  If this variable is not supplied then `dat` will be
#'    treated as a single pedigree.
#' * `indiv`, an individual identifier for each person.  If there are any
#'    duplicated identifiers in the dataset then the family and an underscore
#'    (`_`) will be prepended to all identifiers, and if any duplicates remain
#'    after this then the function will stop executing, with an error message.
#' * `mother`, the individual identifier of each person's mother, or missing
#'    (`NA`) for founders.
#' * `father`, the individual identifier of each person's father, or missing
#'    (`NA`) for founders.
#'
#' @param geno_freq  A vector of strictly positive numbers that sum to `1`.
#' If the possible genotypes of the underlying genetic model are
#' `1:length(geno_freq)` then `geno_freq[j]` is interpreted as the population
#' frequency of genotype `j`, so `geno_freq` is essentially the
#' function `Prior` in the pedigree likelihood on page 117 of (Lange, 2002).
#' For certain genetic models that often occur in applications, these genotype
#' frequencies can be calculated by \code{\link{geno_freq_monogenic}},
#' \code{\link{geno_freq_phased}}, etc.
#'
#' @param trans  An `ngeno^2` by `ngeno` matrix of non-negative numbers whose rows
#' all sum to `1`, where `ngeno = length(geno_freq)` is the number of possible
#' genotypes. The rows of `trans` correspond to joint parental genotypes and
#' the columns correspond to offspring genotypes.  If the possible genotypes
#' are `1:length(geno_freq)` then the element
#' `trans[ngeno * gm + gf - ngeno, go]` is interpreted as the conditional
#' probability that a person has genotype `go`, given that his or her
#' biological mother and father have genotypes `gm` and `gf`, respectively.
#' So `trans` is essentially the transmission function `Tran` on page 117 of
#' (Lange, 2002).  For certain genetic models that often occur in applications,
#' this transmission matrix can be calculated by \code{\link{trans_monogenic}},
#' \code{\link{trans_phased}}, etc.
#'
#' @param penet An `nrow(dat)` by `length(geno_freq)` matrix of non-negative
#' numbers. The element `penet[i,j]` is interpreted as the conditional
#' probability (or probability density) of the phenotype of the person
#' corresponding to row `i` of `dat`, given that his or her genotype is `j`
#' (where the possible genotypes are `1:length(geno_freq)`).
#' Therefore, `penet` is essentially the penetrance function `Pen` on page 117
#' of (Lange, 2002).  If any row of `penet` consists entirely of zeroes then
#' the likelihood is `0`, so the returned log-likelihood will be `-Inf`.
#' Note that genotype data can be incorporated into `penet` by regarding
#' observed genotypes as part of the phenotype, i.e. by regarding observed
#' genotypes as (possibly noisy) measurements of the underlying true genotypes.
#' For example, if the observed genotype of person `i` is `1`
#' (and if genotype measurement error is negligible) then `penet[i,j]`
#' should be `0` for `j != 1` and `penet[i,1]` should be the same as if
#' person `i` were ungenotyped.
#'
#' @param monozyg  An optional list that can be used to specify genetically
#' identical persons, such as monozygotic twins, monozygotic triplets,
#' a monozygotic pair within a set of dizygotic triplets, etc.
#' Each element of the list should be a vector containing the individual
#' identifiers of a group of genetically identical persons, e.g. if `dat`
#' contains six sets of monozygotic twins and one set of monozygotic triplets
#' then `monozyg` will be a list with seven elements, one element a vector of length
#' three and the other six elements all vectors of length two. The order of the list and
#' the orders within its elements do not affect the output of the function.
#' Each group of genetically identical persons should contain two or more
#' persons, the groups should not overlap, and all persons in each group must
#' have the same (non-missing) parents.
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
#' @param load_balancing A logical flag.  When `ncores > 1`, parallelization is
#' achieved either with the function `parallel::parLapply` (if `load_balancing`
#' is `FALSE`) or with the load-balancing function `parallel::parLapplyLB`
#' (if `load_balancing` is `TRUE`, the default). The load-balancing version
#' will usually, but not always, be faster.
#'
#' @details
#' This function provides a fast and general implementation of the
#' Elston-Stewart algorithm to calculate the log-likelihoods of potentially
#' large and complex pedigrees.  General references for the Elston-Stewart
#' algorithm are (Elston & Stewart, 1971), (Lange & Elston, 1975) and
#' (Cannings et al., 1978).
#'
#' Each family within `dat` should be a complete pedigree, meaning that each
#' person should either have both parental identifiers missing (if a founder)
#' or both non-missing (if a non-founder), and each (non-missing) mother or
#' father should have a corresponding row of `dat`.
#'
#' Observed genotypes should be incorporated into `penet`, as described above.
#'
#' The function can handle pedigree loops, such as those
#' caused by inbreeding or by two sisters having children with two brothers
#' from an unrelated family (see (Totir et al., 2009) for a precise definition),
#' though pedigrees with more than a few loops could greatly reduce the speed of
#' the calculation.
#'
#' In `geno_freq`, `trans` and `penet`, the order of the possible genotypes
#' must match, in the sense that the genotype that corresponds to element `j`
#' of `geno_freq` must also correspond to column `j` of `trans` and `penet`,
#' for each `j` in `1:length(geno_freq)`.
#'
#' Sex-specific genetics, such as X-linked genes or genetic loci with sex-specific
#' recombination fractions, can be modelled by letting genotypes `1:nm` be
#' the possible male genotypes and letting `(nm+1):(nm+nf)` be the possible
#' female genotypes, where `nm` and `nf` are the number of possible genotypes
#' for males and females, respectively.  Then, for example, `penet[i,j]` will
#' be `0` if `j %in% 1:nm` and row `i` of `dat` corresponds to a female, and
#' `penet[i,j]` will be `0` if `j %in% (nm+1):(nm+nf)` and row `i` of
#' `dat` corresponds to a male.
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
#' # Settings for a single biallelic locus in Hardy-Weinberg equilibrium
#' # and with a minor allele frequency of 10%
#' geno_freq <- geno_freq_monogenic(c(0.9, 0.1))
#' trans <- trans_monogenic(2)
#'
#' # In dat_small, ora024 and ora027 are identical twins, and so are aey063 and aey064
#' monozyg_small <- list(c("ora024", "ora027"), c("aey063", "aey064"))
#'
#' # Calculate the log-likelihoods for 10 families, each with approximately
#' # 100 family members
#' pedigree_loglikelihood(
#'   dat_small, geno_freq, trans, penet_small, monozyg_small, sum_loglik = FALSE, ncores = 2
#' )
#'
#' # Calculate the log-likelihood for one family with approximately 10,000 family members
#' # Note:  this calculation should take less than a minute on a standard desktop computer
#' # Note:  parallelization would achieve nothing here because there is only one family
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
pedigree_loglikelihood <- function(dat, geno_freq, trans, penet, monozyg = NULL,
                                   sum_loglik = TRUE, ncores = 1, load_balancing = TRUE) {

  #check_dat(dat)

  # Check if family can be found
  one_fam <- FALSE
  if (!"family" %in% names(dat)) {
    one_fam <- TRUE
    dat <- cbind(family = 1, dat)
  }

  # Deal with identical twins, etc
  if (length(monozyg) >= 1) {
    for (i in 1:length(monozyg)) {
      birth_group <- which(dat$indiv %in% monozyg[[i]])
      if (length(birth_group) >= 2) {
        representative <- birth_group[1]
        p <- apply(penet[birth_group,], 2, prod)
        penet[representative,] <- p
        dat$mother[!is.na(dat$mother) & dat$mother %in% monozyg[[i]]] <- dat$indiv[representative]
        dat$father[!is.na(dat$father) & dat$father %in% monozyg[[i]]] <- dat$indiv[representative]
        birth_group <- setdiff(birth_group, representative)
        penet <- penet[-birth_group,]
        dat <- dat[-birth_group,]
      }
    }
  }

  # Check if IDs are unique, if not then concatenate with family ID
  if (one_fam & any(duplicated(dat$indiv)))  stop("Individual identifiers are not unique")
  if (!one_fam & any(duplicated(dat$indiv))) {
    ID <- paste0(dat$family, "_", dat$indiv)
    if (any(duplicated(ID)))  stop("Individual identifiers cannot be made unique")
    mumID <- paste0(dat$family, "_", dat$mother)
    mumID[is.na(dat$mother)] <- NA
    dadID <- paste0(dat$family, "_", dat$father)
    dadID[is.na(dat$father)] <- NA
    dat$indiv <- ID
    dat$mother <- mumID
    dat$father <- dadID
  }

  # Convert IDs to consecutive integers
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
    if (length(keep) == 1)  peneti <- matrix(peneti, 1, )
    fam_penet_list[[i]] <- list(
      fam = fam,
      penet = peneti
    )
  }
  names(fam_penet_list) <- ufamID
  #fam_penet <- fam_penet_list[[1]]

  if (ncores == 1) {
    ll <- lapply(fam_penet_list, pedigree_loglikelihood_g,
                 geno_freq = geno_freq, trans = trans)
  } else {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (load_balancing) {
      ll <- parallel::parLapplyLB(
        cl = cl, X = fam_penet_list, fun = pedigree_loglikelihood_g,
        geno_freq = geno_freq, trans = trans)
    } else {
      ll <- parallel::parLapply(
        cl = cl, X = fam_penet_list, fun = pedigree_loglikelihood_g,
        geno_freq = geno_freq, trans = trans)
    }
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
