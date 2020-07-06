#' Simulated data on one family with approximately 10,000 members.
#'
#' A dataset giving the relationship structure of one large family
#' and phenotypic data on the family members
#'
#' @format A data frame with 10,002 rows (corresponding to persons) and 6 variables:
#' \describe{
#'   \item{indiv}{an identifier (ID) for the person}
#'   \item{mother}{the individual ID of the person's mother}
#'   \item{father}{the individual ID of the person's father}
#'   \item{sex}{the person's sex (1 = male, 2 = female)}
#'   \item{aff}{the person's disease status (1 = case, 0 = control)}
#'   \item{age}{the person's age, in years}
#' }
#' @source Simulated
"dat_large"
