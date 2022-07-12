## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
library(clipp)
MAF <- 0.1
geno_freq <- geno_freq_monogenic(p_alleles = c(1 - MAF, MAF))
trans <- trans_monogenic(n_alleles = 2)

## -----------------------------------------------------------------------------
geno_freq_monogenic(p_alleles = c(1 - MAF, MAF), annotate = TRUE)

## -----------------------------------------------------------------------------
trans_monogenic(n_alleles = 2, annotate = TRUE)

## -----------------------------------------------------------------------------
data("dat_small", "penet_small", "dat_large", "penet_large")
head(dat_small)
head(penet_small)

## -----------------------------------------------------------------------------
monozyg_small <- list(c("ora024", "ora027"), c("aey063", "aey064"))

## -----------------------------------------------------------------------------
pedigree_loglikelihood(dat_small, geno_freq, trans, penet_small, 
                       monozyg = monozyg_small, sum_loglik = FALSE, ncores = 2)

## ---- eval = FALSE------------------------------------------------------------
#  system.time(ll <- pedigree_loglikelihood(dat_large, geno_freq, trans, penet_large))
#  #>    user  system elapsed
#  #>   10.64    0.15   10.83
#  ll
#  #> [1] -18020.99

## -----------------------------------------------------------------------------
head(dat_small)

## -----------------------------------------------------------------------------
penet.fn <- function(i, logodds1, logodds2) {
  prob1 <- 1/(1 + exp(-logodds1))
  prob2 <- 1/(1 + exp(-logodds2))
  penet.i <- c(prob1, prob1, prob2)
  if (dat_small$aff[i] == 0)  penet.i <- 1 - penet.i
  return(penet.i)
}

## ---- eval = FALSE------------------------------------------------------------
#  library(stats4)
#  minusll <- function(logodds1 = 0, logodds2 = 0) {
#    penet <- t(sapply(1:nrow(dat_small), penet.fn, logodds1, logodds2))
#    loglik <- pedigree_loglikelihood(dat_small, geno_freq, trans, penet,
#                                     monozyg = monozyg_small, ncores = 2)
#    return(-loglik)
#  }
#  minusll()
#  #> [1] 705.6238
#  fit <- mle(minusll)
#  summary(fit)
#  #> Maximum likelihood estimation
#  #>
#  #> Call:
#  #> mle(minuslogl = minusll)
#  #>
#  #> Coefficients:
#  #>           Estimate Std. Error
#  #> logodds1 -1.359962  0.1054336
#  #> logodds2 -1.313467  6.8597364
#  #>
#  #> -2 log L: 1030.901

## -----------------------------------------------------------------------------
head(dat_small)

## -----------------------------------------------------------------------------
penet.fn <- function(i, logodds1, logodds2) {
  prob1 <- 1/(1 + exp(-logodds1))
  prob2 <- 1/(1 + exp(-logodds2))
  penet.i <- c(prob1, prob1, prob2)
  if (dat_small$aff[i] == 0)  penet.i <- 1 - penet.i
  if (dat_small$geno[i] == "1/1")  penet.i[-1] <- 0     ###
  if (dat_small$geno[i] == "1/2")  penet.i[-2] <- 0     ###
  if (dat_small$geno[i] == "2/2")  penet.i[-3] <- 0     ###
  return(penet.i)
}

## ---- eval = FALSE------------------------------------------------------------
#  minusll()
#  #> [1] 788.5003
#  fit <- mle(minusll)
#  summary(fit)
#  #> Maximum likelihood estimation
#  #>
#  #> Call:
#  #> mle(minuslogl = minusll)
#  #>
#  #> Coefficients:
#  #>           Estimate Std. Error
#  #> logodds1 -1.357596 0.08405912
#  #> logodds2 -1.395321 0.61632248
#  #>
#  #> -2 log L: 1196.65
#  

## -----------------------------------------------------------------------------
head(dat_small)
genotype_probabilities(target = "ora008", dat_small, geno_freq, trans, 
                       penet_small, monozyg_small)

## -----------------------------------------------------------------------------
penet.fn <- function(i) {
  penet.i <- rep(1, 3)
  if (dat_small$geno[i] == "1/1")  penet.i[-1] <- 0
  if (dat_small$geno[i] == "1/2")  penet.i[-2] <- 0
  if (dat_small$geno[i] == "2/2")  penet.i[-3] <- 0
  return(penet.i)
}
penet <- t(sapply(1:nrow(dat_small), penet.fn))
genotype_probabilities(target = "ora008", dat_small, geno_freq, trans, 
                       penet, monozyg_small)

