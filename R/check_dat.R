check_dat <- function(dat) {

  # Check if the input is a data.frame
  if (class(dat) != "data.frame") {
    #stop(paste0(deparse(substitute(dat)), " is not a data.frame"))
    stop(sprintf("%s is not a data.frame", deparse(substitute(dat))),
         call. = FALSE)
  }

  # Check if all required columns are here
  column_names <- c("indiv", "mother", "father")
  w <- which(!column_names %in% names(dat))
  if (!all(column_names %in% names(dat))) {
    stop(paste("Cannot find column", paste(sQuote(column_names[w]),
                                            collapse = ", ")), call. = FALSE)
  }

  # Get columns
  indiv <- dat$indiv
  mother <- dat$mother
  father <- dat$father

  # Check if any missing values in column 'indiv'
  if (anyNA(indiv)) {
    stop("Column 'indiv' should have no missing values", call. = FALSE)
  }

  # Check if 'indiv' column has unique values
  if (sum(duplicated(indiv)) >= 1) {
    stop("Column 'indiv' should be unique for each person", call. = FALSE)
  }

}

