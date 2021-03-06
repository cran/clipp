% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dat_large.R
\docType{data}
\name{dat_large}
\alias{dat_large}
\title{Simulated data on one family with approximately 10,000 members}
\format{
A data frame with 10,002 rows (corresponding to persons) and the
following 6 variables:
\describe{
\item{\code{indiv}}{an individual identifier for each person}
\item{\code{mother}}{the individual identifier of the person's mother}
\item{\code{father}}{the individual identifier of the person's father}
\item{\code{sex}}{the person's sex (1 = male, 2 = female)}
\item{\code{aff}}{the person's disease status (1 = case, 0 = control)}
\item{\code{age}}{the person's last known age in years (for controls) or
age at diagnosis in years (for cases)}
}
}
\source{
Simulated
}
\usage{
dat_large
}
\description{
A dataset giving the relationship structure of one large family
and phenotypic data on the family members
}
\keyword{datasets}
