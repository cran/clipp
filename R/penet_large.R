#' A penetrance matrix relating the phenotypes in `dat_large` to three genotypes
#'
#' A matrix relating the phenotypes of `dat_large` to the three unphased
#' genotypes of a single biallelic, autosomal genetic locus.
#' The element `penet_large[i,j]` is the conditional probability of
#' the phenotypes (i.e. `sex`, `aff` and `age`) of the person
#' in row `i` of `dat_large`, given that
#' his or her genotype is `j` (here labelling the genotypes as `1, 2, 3`,
#' where genotype `2` is the heterozygous genotype).
#'
#' @format A matrix with 10,002 rows (corresponding to persons) and 3
#' columns (corresponding to genotypes).
#'
#' @source Simulated
"penet_large"
