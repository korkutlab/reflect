# z-score transformation -------
#' @export
zscore <- function (vec) {
  vec <- (vec - mean(vec, na.rm = TRUE)) / sd(vec, na.rm = TRUE)
  return(vec)
}

# median center transformation -------
#' @export
median_center <- function(vec) {
  vec <- vec - median(vec)
  return(vec)
}