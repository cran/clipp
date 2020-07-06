#' Simulated data on 10 families with approximately 100 members each.  
#' 
#' A dataset giving the relationship structure of 10 families 
#' and phenotypic data on the family members
#'
#' @format A data frame with 1018 rows (corresponding to persons) and 7 variables:
#' \describe{
#'   \item{family}{an identifier (ID) for the person's family}
#'   \item{indiv}{an ID for the person}
#'   \item{mother}{the individual ID of the person's mother}
#'   \item{father}{the individual ID of the person's father}
#'   \item{sex}{the person's sex (1 = male, 2 = female)}
#'   \item{aff}{the person's disease status (1 = case, 0 = control)}
#'   \item{age}{the person's age, in years}
#' }
#' @source Simulated
"dat_small"