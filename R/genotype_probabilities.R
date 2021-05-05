#' Calculate genotype probabilities for a target person
#'
#' For a chosen individual within a specified family, calculate the person's
#' conditional genotype probabilities, given the family's phenotypes and
#' relationship structure
#'
#' @param target  The individual identifier (an element of `fam$indiv`)
#' of the person in the pedigree `fam` whose genotype probabilities are being
#' sought.
#'
#' @param fam A data frame specifying the family's relationship structure,
#' with rows corresponding to people and columns
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
#' If the possible genotypes of the underlying genetic model are
#' `1:length(geno_freq)` then `geno_freq[j]` is interpreted as the population
#' frequency of genotype `j`.
#' For certain genetic models that often occur in applications, these genotype
#' frequencies can be calculated by \code{\link{geno_freq_monogenic}},
#' \code{\link{geno_freq_phased}}, etc.
#'
#' @param trans  An `ngeno^2` by `ngeno` matrix of non-negative numbers whose rows
#' all sum to `1`, where `ngeno = length(geno_freq)` is the number of possible
#' genotypes.  The rows of `trans` correspond to joint parental
#' genotypes and the columns correspond to offspring genotypes.
#' If the possible genotypes are `1:length(geno_freq)` then the element
#' `trans[ngeno * gm + gf - ngeno, go]` is interpreted as the conditional
#' probability that a person has genotype `go`, given that his or her
#' biological mother and father have genotypes `gm` and `gf`, respectively.
#' For certain genetic models that often occur in applications,
#' this transmission matrix can be calculated by \code{\link{trans_monogenic}},
#' \code{\link{trans_phased}}, etc.
#'
#' @param penet An `nrow(fam)` by `length(geno_freq)` matrix of non-negative
#' numbers. The element `penet[i,j]` is interpreted as the conditional
#' probability (or probability density) of the phenotype of the person
#' corresponding to row `i` of `fam`, given that his or her genotype is `j`
#' (where the possible genotypes are `1:length(geno_freq)`).
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
#' identifiers of a group of genetically identical persons, e.g. if `fam`
#' contains a set of monozygotic twins (and no other genetically identical
#' persons) then `monozyg` will be a list
#' with one element, and that element will be a vector of length two containing
#' the individual identifiers of the twins.  The order of the list
#' and the orders of its elements do not affect the output of the function.
#' Each group of genetically identical persons should contain two or more
#' persons, the groups should not overlap, and all persons in each group must
#' have the same (non-missing) parents.
#'
#' @details
#' The genotype probabilities are calculated by essentially the same algorithm
#' as \code{\link{pedigree_loglikelihood}};
#' see there for details.  The genotype probabilities only depend on the
#' connected component of the pedigree that contains `target`, so the
#' function first restricts `fam` and `penet` to the rows corresponding to this
#' connected component.  For example, if `fam` is the union of two
#' unrelated families then this function will restrict to the subfamily
#' containing `target` before performing the calculation.
#'
#' @return A vector of length `length(geno_freq)`, whose `j`th element is
#' the conditional probability that the target person has genotype `j`,
#' given the family's relationship structure and phenotypes.  A vector of `NA`s
#' will be returned if a row of `penet` consists entirely of zeroes or if
#' the pedigree is impossible for any other reason
#' (after restricting `fam` and `penet` to the connected component of
#' the pedigree containing `target`).
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
#' monozyg <- list(c("ora024", "ora027"))  # ora024 and ora027 are identical twins
#' trans <- trans_monogenic(2)
#' geno_freq <- geno_freq_monogenic(p_alleles = c(0.9, 0.1))
#' genotype_probabilities(target = "ora008", fam, geno_freq, trans, penet, monozyg)
#'
genotype_probabilities <- function(target, fam, geno_freq, trans, penet,
                                   monozyg = NULL) {

  dat <- fam[, names(fam) != "family"]
  dat <- data.frame(family = rep(1,nrow(fam)), dat)

  if (is.null(target) | !target %in% dat$indiv) {
    stop("The target cannot be found in the family data",  call. = FALSE)
  }

  # Deal with the case where the target forms a singleton subfamily of dat
  parents <- setdiff(c(dat$mother, dat$father), NA)
  keep <- dat$indiv == target
  if (is.na(dat$mother[keep]) & is.na(dat$father[keep]) & !(target %in% parents)) {
    p <- geno_freq * penet[keep,]
    if (sum(p) == 0) {return(rep(NA,length(geno_freq)))} else {return(p/sum(p))}
  }

  # Restrict the family to its connected component containing target
  # (without this, whether or not you get NAs as output is slightly random)
  nold <- 0
  keepIDs <- target
  while (nold < length(keepIDs)) {
    nold <- length(keepIDs)
    new <- dat$indiv[dat$mother %in% keepIDs]
    new <- c(new, dat$indiv[dat$father %in% keepIDs])
    new <- c(new, dat$mother[dat$indiv %in% keepIDs])
    new <- c(new, dat$father[dat$indiv %in% keepIDs])
    new <- setdiff(new, NA)
    keepIDs <- unique(c(keepIDs, new))
  }
  keep <- dat$indiv %in% keepIDs
  dat <- dat[keep,]
  penet <- penet[keep,]

  # Deal with identical twins, etc
  if (length(monozyg) >= 1) {
    for (i in 1:length(monozyg)) {
      birth_group <- which(dat$indiv %in% monozyg[[i]])
      if (length(birth_group) >= 2) {
        if (target %in% monozyg[[i]]) {
          representative <- which(dat$indiv==target)
        } else {
          representative <- birth_group[1]
        }
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

  # Save which row corresponds to target because IDs will be converted to
  # row numbers soon
  target_id <- which(dat$indiv == target)

  # Convert IDs
  dat$mother[which(is.na(dat$mother))] <- ""
  dat$father[which(is.na(dat$father))] <- ""
  dat <- convert_IDs(dat, convert.IDs.numeric = TRUE)

  # Perform the calculation and output
  fam_penet <- list(fam = dat, penet = penet)
  p <- pedigree_loglikelihood_g(fam_penet, geno_freq, trans, target_id)
  return(p)

}

