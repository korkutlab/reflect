#' @export
print.refdata <- function(x,...) {
  cat("reflect data object", fill = TRUE)
  cat(fill = TRUE)

  cat(paste0(names(x)[1], ": matrix for expression/alteration values"), fill = TRUE)
  cat("\t", "Row: SampleID",
      paste0("(#samples = ", nrow(x[[1]]), ")"),
      fill = TRUE)
  cat("\t", "Col: Feature",
      paste0("(#features = ", ncol(x[[1]]), ")"),
      fill = TRUE)

  cat(paste0(names(x)[2], ": data frame for sample annotations, having"), ncol(x[[2]]), "columns:", fill = TRUE)
  for (name in colnames(x[[2]])) {
    cat("\t", name, fill = TRUE)
  }

  cat(paste0(names(x)[3], ": data frame for feature annotations, having"), ncol(x[[3]]), "columns:", fill = TRUE)
  for (name in colnames(x[[3]])) {
    cat("\t", name, fill = TRUE)
  }
}


#' @export
print.refresult <- function(x,...) {
  cat("reflect result object", fill=TRUE)
  cat(fill=TRUE)

  cat(paste0(names(x)[1], ": Gap statistic and best wbound."), fill = TRUE)
  cat(paste0(names(x)[2], ": An object of class hclust."), fill = TRUE)
  cat(paste0(names(x)[3], ": Recurrence P values."), fill = TRUE)
  cat(paste0(names(x)[3], ": Rcurrent and acitionable features."), fill = TRUE)
  cat(paste0(names(x)[5], ": Co-altered, recurrent, actionable targets."), fill = TRUE)
}


#' @export
print.refcolor <- function(x,...) {
  cat("reflect color object", fill=TRUE)
  cat(fill=TRUE)

  cat(paste0(names(x)[1], ": Colors for TCGA tumor types."), fill = TRUE)
  cat(paste0(names(x)[2], ": Colors for cell line tumor types."), fill = TRUE)
  cat(paste0(names(x)[3], ": Colors for categories."), fill = TRUE)
}

