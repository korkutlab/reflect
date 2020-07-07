#' REFLECT
#'
#' REcurrent Features LEveraged for Combination Therapies
#'
#' REFLECT is to accelerate discovery of combination therapies tailored to multi-omic profiles of cancer patient cohorts with matching oncogenic features. REFLECT identifies co-actionable, co-occurring oncogenic alterations that are recurrent within patient cohorts using genomic, transcriptomic and phosphoproteomic data. The REFLECT predictions are testable hypotheses on responses to drug combinations in specific patient cohorts matched to co-actionable alterations. The analysis of proteomic and transcriptomic signatures in relation to the variant information in patient cohorts also enable prediction of functional rare variants. REFLECT provides a novel, robust and widely applicable precision medicine approach and resource for matching oncogenic alterations to combination therapies.
#'
#' @name reflect-package
#' @aliases reflect-package reflect
#' @docType package
#' @references
#' Li X., et al. (2020) Precision combination therapies from recurrent oncogenic co-alterations \doi{https://doi.org/10.1101/2020.06.03.132514}.
#' @keywords package
NULL


#' Dataset of a EGFR-mutant Cohort
#'
#' This dataset is (phospho-)protein expression Z scores of an EGFR-mutant cohort measured by RPPA.
#'
#' @name egfr_data
#' @docType data
#' @format The format is a list containing the following elements:
#' - mat_value: (Phospho-)protein expression Z scores from RPPA mesurement (260 samples, 179 features).
#' - df_sample: A data frame for sample annotations, having 4 columns: SampleID, TumorType, Stratification, Category.
#' - df_feature: A data frame for feature annotations, having 3 columns: Feature, FunctionScore, IsActionable.
#'
#' @references
#'
#' Li X., et al. (2020) Precision combination therapies from recurrent oncogenic co-alterations \doi{https://doi.org/10.1101/2020.06.03.132514}.
#'
#' @keywords datasets
#' @examples
#' egfr_data
#'
NULL


#' Result of a EGFR-mutant Cohort from REFLECT
#'
#' This is the result from running reflect_pipeline on the egfr_data dataset.
#'
#' @name egfr_result
#' @docType data
#' @format The format is a list containing the following elements:
#' - gapstat_bestwbound: Gap statistic and best wbound.
#' - shc: An object of class hclust.
#' - mat_recur_pval: A matrix of recurrence P values.
#' - mat_recur_pval: Rcurrent and acitionable features.
#' - df_coaltered_targets: Co-altered, recurrent, actionable targets.
#'
#' @references
#'
#' Li X., et al. (2020) Precision combination therapies from recurrent oncogenic co-alterations \doi{https://doi.org/10.1101/2020.06.03.132514}.
#'
#' @keywords datasets
#' @examples
#' egfr_result
#'
NULL


#' Color for Categories from REFLECT
#'
#' This is a collection of colors for tumor types and variant categories.
#'
#' @name reflect_color
#' @docType data
#' @format The format is a list containing the following elements:
#' - df_tcga_color: Colors for TCGA tumor types.
#' - df_ccl_color: Colors for cell line tumor types.
#' - df_category_color: Colors for categories.
#'
#' @references
#'
#' Li X., et al. (2020) Precision combination therapies from recurrent oncogenic co-alterations \doi{https://doi.org/10.1101/2020.06.03.132514}.
#'
#' @keywords datasets
#' @examples
#' reflect_color
#'
NULL


#' Pipe
#'
#' Pipe an object forward into a function or call expression.
#'
#' @importFrom dplyr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs	A value or the magrittr placeholder.
#' @param rhs	A function call using the magrittr semantics.
#' @examples
#' c(5:1) %>% sort()
#'
NULL

