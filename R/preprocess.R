# z-score transformation -------
#' Standardize a numerical vector (Z score transformation)
#'
#' Given a numerical vector x, return (x - mean(x)) / sd(x), where NA are removed in the calculation of mean(x) and sd(x).
#'
#' @param x A numberical vector (at least 2 non-NA items).
#' @return The standardized form of \code{x}.
#' @examples
#' x <- 1:10
#' zscore(x)
#' @export
zscore <- function (x) {
  if (sum(!is.na(x)) < 2) stop("x must have at least 2 non-NA items")
  x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  return(x)
}

# median center transformation -------
#' Put the median of a numerical vector to 0
#'
#' Given a numerical vector x, return x - median(x), where NA are removed in the calculation of median(x)
#'
#' @param x A numberical vector (at least 2 non-NA items).
#' @return The zero-median form of \code{x}.
#' @examples
#' x <- 1:10
#' median_center(x)
#' @export
median_center <- function(x) {
  if (sum(!is.na(x)) < 2) stop("x must have at least 2 non-NA items")
  x <- x - median(x)
  return(x)
}
