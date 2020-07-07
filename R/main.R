# best_wbound function -------
#' Get the Best Tunning Parameter
#'
#' This employs a permutation approach to select the tuning parameter, which controls sparsity of features through L1 regularization. Note this calculation may take a long time.
#'
#' @param mat_value A matrix of expression/alteration with samples as rows and features as columns.
#' @param dissimilarity A string for the type of dissimilarity, either "squared.distance" or "absolute.value". Default "squared.distance".
#' @param wbounds The sequence of tuning parameters to consider. If NULL, then a default sequence seq(1.1, sqrt(ncol(mat_value)), 100) will be used. If non-null, should be greater than 1.
#' @param nperms The number of permutations to perform. Default 10.
#' @param min_number_features The minimal number of features that the best wbound could generate. Only wbounds that generates more than this number of features can considered. Default 10.
#' @return \item{df_gapstat}{A data frame that contains columns: Wbound, NumberFeatures, GapStat, SD.} \item{best_wbound}{The best tunning parameter whose gap statistic is within one-standard-error to the left of the maximum gap statistic.}
#' @examples
#' library(reflect)
#' mat_value <- egfr_data$mat_value
#' gapstat_bestwbound <- get_best_wbound(mat_value)
#' gapstat_bestwbound$best_wbound
#' @export
get_best_wbound <- function(mat_value,
                            dissimilarity = c("squared.distance",
                                              "absolute.value"),
                            wbounds = NULL,
                            nperms = 10,
                            min_number_features = 10
) {
  if (sum(duplicated(rownames(mat_value))) > 0 |
      sum(duplicated(colnames(mat_value))) > 0) {
    stop('samples and features (rownames and colnames of mat_value) need to be unique')
  }

  dissimilarity <- match.arg(dissimilarity)

  if (is.null(wbounds)) {
    wbounds <- seq(1.1, sqrt(ncol(mat_value)) * 1.0, len = 100)
  }

  perm.out <- sparcl::HierarchicalSparseCluster.permute(mat_value,
                                                        dissimilarity = "squared.distance",
                                                        wbounds = wbounds,
                                                        nperms = nperms)
  tb <- dplyr::tibble(perm.out$wbounds, perm.out$nnonzerows, perm.out$gaps, perm.out$sdgaps)
  colnames(tb) <- c("Wbound", "NumberFeatures", "GapStat", "SD")
  df_gapstat <- as.data.frame(tb)

  # lower bound contraint
  tb <- tb %>%
    dplyr::filter(NumberFeatures >= min_number_features)

  # keep maximum gap statistic for each number of features
  tb <- tb %>%
    dplyr::arrange(NumberFeatures, desc(GapStat)) %>%
    dplyr::distinct(NumberFeatures, .keep_all = TRUE)
  tb_gap <- tb

  # get best_wbound ------
  gap_max <- max(tb_gap$GapStat)
  tb <- tb_gap[c(1:which.max(tb_gap$GapStat)), ]
  gap_mean <- mean(tb$GapStat)
  gap_sd <- sd(tb$GapStat)
  gap_se <- gap_sd/sqrt(nrow(tb))
  while (nrow(tb) > 1) {
    gap_below <- tb$GapStat < gap_mean - gap_sd
    if (sum(gap_below) == 0) gap_below <- tb$GapStat < gap_mean
    beg <- max(which(gap_below)) + 1
    tb <- tb[c(beg:nrow(tb)), ]  # only keep a coutinuous region near the maximum
    if (nrow(tb) < 2) break()
    gap_mean_old <- gap_mean
    gap_sd_old <- gap_sd
    gap_se_old <- gap_se
    gap_mean <- mean(tb$GapStat)
    gap_sd <- sd(tb$GapStat)
    gap_se <- gap_sd/sqrt(nrow(tb))
    if (abs(gap_mean - gap_mean_old) < gap_se_old) break()
  }
  best_wbound <- tb$Wbound[1]

  return(list("df_gapstat" = df_gapstat,
              "best_wbound" = best_wbound))
}


# sparse_hclust -----
#' Run a Sparse Hierarchical Clustering
#'
#' Given a matrix of expression/alteration and the wbound value for L1 regularization, it performs a sparse hierarchical clustering.
#'
#' @param mat_value A matrix of expression/alteration with samples as rows and features as columns.
#' @param wbound A real number as the wbound used in sparse hierarchial clustering.
#' @param dissimilarity A string for the type of dissimilarity, either "squared.distance" or "absolute.value". Default "squared.distance".
#' @return \item{hc}{An object of class hclust which describes the tree produced by the clustering process. See detail in function hclust from stats.} \item{weight}{The weights of features used in sparse hierachical clustering.} \item{features_nonzero}{Features that have nonzero weights.} \item{mat_value_clustered}{A matrix of values that are clustered by the sparse hierarchical clustering.}
#' @examples
#' library(reflect)
#' mat_value <- egfr_data$mat_value
#' wbound <- 2.0
#' shc <- sparse_hclust(mat_value, wbound)
#' @export
sparse_hclust <- function(mat_value,
                          wbound,
                          dissimilarity = c("squared.distance",
                                            "absolute.value")
) {
  if (sum(duplicated(rownames(mat_value))) > 0 |
      sum(duplicated(colnames(mat_value))) > 0) {
    stop('samples and features (rownames and colnames of mat_value) need to be unique')
  }

  dissimilarity <- match.arg(dissimilarity)

  sparsehc <- sparcl::HierarchicalSparseCluster(x = mat_value,
                                                dissimilarity = dissimilarity,
                                                dists = NULL,
                                                wbound = wbound)

  # the original hc has numbers as lables
  # change the labels to sampleid
  sparsehc$hc$labels <- rownames(mat_value)

  weight <- sparsehc$ws
  names(weight) <- colnames(mat_value)
  weight <- sort(weight, decreasing = TRUE)

  features_nonzero <- names(weight)[weight > 0]
  mat_value_clustered <- mat_value[sparsehc$hc$order, features_nonzero]

  return(list("hc" = sparsehc$hc,
              "weight" = weight,
              "features_nonzero" = features_nonzero,
              "mat_value_clustered" = mat_value_clustered))
}


# robustness test -------
get_top_n_features <- function(mat_value, n_top,
                               dissimilarity, wbounds, nperms,
                               min_number_features
) {
  gapstat_bestwbound <- get_best_wbound(mat_value,
                                   dissimilarity = dissimilarity,
                                   wbounds = wbounds,
                                   nperms = nperms,
                                   min_number_features = min_number_features)
  best_wbound <- gapstat_bestwbound$best_wbound
  sparsehc <- sparcl::HierarchicalSparseCluster(x = mat_value,
                                                dissimilarity = dissimilarity,
                                                dists = NULL,
                                                wbound = best_wbound)
  tb <- dplyr::tibble(Feature = colnames(mat_value), ws = sparsehc$ws[, 1])
  tb <- tb %>%
    dplyr::filter(ws > 0) %>%
    dplyr::arrange(dplyr::desc(ws))
  if(nrow(tb) > n_top) tb <- tb[1:n_top, ]
  return(tb$Feature)
}

get_score <- function(mat_value, indices, features_ref, n_top_sub,
                      dissimilarity, wbounds, nperms,
                      min_number_features){
  mat_sub <- mat_value[indices, ]
  features_sub <- get_top_n_features(mat_sub, n_top_sub,
                                     dissimilarity, wbounds, nperms,
                                     min_number_features)
  score <- 0
  for (f in features_sub) {
    pos <- match(f, features_ref)
    if (!is.na(pos)) {
      if (pos <= n_top_sub) {
        score <- score + 1
      } else if (pos <= 2*n_top_sub) {
        score <- score + (2*n_top_sub - pos)/n_top_sub
      }
    }
  }
  return(score)
}

#' Robustness Analysis Due to Subsampling
#'
#' This gives a measure of how robust the sparse hierarchical clustering performs on a cohort by subsampling the samples from the target cohort.
#'
#' @param mat_value A matrix of expression/alteration with samples as rows and features as columns.
#' @param dissimilarity A string for the type of dissimilarity, either "squared.distance" or "absolute.value". Default "squared.distance".
#' @param wbounds The sequence of tuning parameters to consider. If NULL, then a default sequence seq(1.1, sqrt(ncol(mat_value)), 100) will be used. If non-null, should be greater than 1.
#' @param nperms The number of permutations to perform. Default 10.
#' @param min_number_features The minimal number of features that the best wbound could generate. Only wbounds that generates more than this number of features can considered. Default 10.
#' @param n_top_ref The number of top features used as as reference from the full cohort. Default 20.
#' @param n_top_sub The number of top features to be compared from the subsampled cohort. Default 10.
#' @param fs_subsample The sequence of fractions to be subsampled. Default seq(0.6, 0.95, 0.05).
#' @param n_times The number of sampling times for each fraction. Default 200.
#' @param no_cores This function can be run in parallel, in which no_cores is the number of cores. If NULL, (the number of all available cores - 1) is used.
#' @param seed The random seed used in subsampling. Default 123.
#' @return A matrix that has columns: \item{f_subsample}{Fraction of data used.} \item{score}{Concordance between the features from partial data and the reference features from the full data.}
#' @examples
#' library(reflect)
#' mat_value <- egfr_data$mat_value
#' mat_robust <- subsample_robustness(mat_value, no_cores = 1) # this may take a long time if only 1 core is used.
#' @export
subsample_robustness <- function(mat_value,
                                 dissimilarity = c("squared.distance",
                                                   "absolute.value"),
                                 wbounds = NULL,
                                 nperms = 10,
                                 min_number_features = 10,
                                 n_top_ref = 20,
                                 n_top_sub = 10,
                                 fs_subsample = seq(0.6, 0.95, 0.05),
                                 n_times = 200,
                                 no_cores = NULL,
                                 seed = 123
) {
  if (sum(duplicated(rownames(mat_value))) > 0 |
      sum(duplicated(colnames(mat_value))) > 0) {
    stop('samples and features (rownames and colnames of mat_value) need to be unique')
  }

  dissimilarity <- match.arg(dissimilarity)

  # create reference for n_top_ref from the full dataset
  features_ref <- get_top_n_features(mat_value, n_top_ref,
                                     dissimilarity, wbounds, nperms,
                                     min_number_features)

  # subsampling and calculate robustness score ------
  n_totsample <- nrow(mat_value)

  indices_list <- list()
  f_subsample_vec <- c()
  set.seed(seed)
  for (f_subsample in fs_subsample) {
    for (i in c(1:n_times)) {
      indices = sample(1:n_totsample, size=floor(n_totsample*f_subsample))
      f_subsample_vec <- c(f_subsample_vec, f_subsample)
      indices_list <- c(indices_list, list(indices))
    }
  }

  # parallel
  if (is.null(no_cores)) no_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(no_cores)

  parallel::clusterEvalQ(cl, library("tidyverse"))
  parallel::clusterEvalQ(cl, library("sparcl"))
  parallel::clusterExport(cl=cl,
                          varlist=c("get_best_wbound", "get_top_n_features", "get_score",
                                    "indices_list", "mat_value", "features_ref", "n_top_sub",
                                    "dissimilarity", "wbounds", "nperms",
                                    "min_number_features"),
                          envir=environment())

  get_score_wrap <- function(id){
    indices <- indices_list[[id]]
    get_score(mat_value, indices, features_ref, n_top_sub,
              dissimilarity, wbounds, nperms,
              min_number_features)
  }

  scores <- c()
  id_beg <- 1
  while (id_beg <= length(indices_list)) {
    id_end <- min(id_beg+no_cores-1, length(indices_list))
    message('Running subprocesses ', id_beg, ' to ', id_end, ' (total processes: ', length(indices_list), ')')

    s <- parallel::parSapply(cl, c(id_beg:id_end), get_score_wrap)
    scores <- c(scores, s)
    id_beg <- id_beg + no_cores
  }

  parallel::stopCluster(cl)

  # combine results
  tb <- dplyr::tibble(f_subsample = f_subsample_vec,
                      score = scores) %>%
    mutate(score = score/n_top_sub)

  return(as.matrix(tb))
}


# functions of recurrence test ------
# pval calculation
consecutive_hits_pval <- function(p, n, k) {
  q <- 1 - p

  # check p
  if (p == 0) {
    return(NA)
  } else if (p == 1) {
    return(1)
  }

  # check k
  if (k < 1) {
    message("k must be a postive integer\n")
    stop()
  }
  # trival case when k = 1
  if (k == 1) {
    return(1 - q**n)
  }

  # dealing case k > 1

  # initialize
  P <- replicate(n, 0)
  P0 <- replicate(n+1, 0)

  # note P is defined in 1:n; P0 is defined in 2:n+1
  for (t in c(1:(k-1))) {
    P[t] <- 1
  }
  P[k] <- 1 - p**k
  for (t in c(2:(k+1))) {
    P0[t] <- q*P[t-1]
  }

  # since the recurrence of Q0 may introduce numerical sigularity,
  # we use the recurrence of P0
  # iterate P0 from k+2 to n+1
  for (t in c((k+2):(n+1))) {
    for (i in c(0:(k-1))) {
      P0[t] <- P0[t] + p**i*P0[t-1-i]
    }
    P0[t] <- q*P0[t]
  }

  # calculate P
  for (t in c((k+1):n)) {
    P[t] <- P0[t+1]/q
  }
  return(max(1 - P[n], 0)) # it may be negative due to floating point errors, so use 0 as a cutoff
}

# recurrence pval fun
recur_pval <- function(vec) {
  num_sel <- sum(vec > 0)
  num_total <- length(vec)

  vec_pval <- rep(1, length(vec))
  names(vec_pval) <- names(vec)

  idx_box <- c()
  idx <- 1L
  while (idx <= length(vec)) {
    if (vec[idx] > 0) {
      idx_box <- c(idx_box, idx)
    }
    if (!(vec[idx] > 0) | idx == length(vec)) {
      num <- length(idx_box)
      if (num > 0L) {
        # pval calculation
        vec_pval[idx_box] <- consecutive_hits_pval(num_sel/num_total, num_total, num)

        # empty idx box
        idx_box <- c()
      }
    }
    idx <- idx + 1L
  }
  return(vec_pval)
}

# get recurrence pval
get_recur_pval_kernel <- function(mat_value, vec_fs) {
  # score matrix
  mat_fs <- matrix(rep(vec_fs, nrow(mat_value)),
                   ncol = length(vec_fs), byrow = TRUE)
  colnames(mat_fs) <- names(mat_value)
  rownames(mat_fs) <- rownames(mat_value)
  mat_score <- mat_value * mat_fs

  # discreted score; we are interested in positive score
  cutoff_mag <- 0
  mat_score_discrete <- matrix(0L,
                               nrow = nrow(mat_score),
                               ncol = ncol(mat_score),
                               dimnames = dimnames(mat_score))
  mat_score_discrete[mat_score > cutoff_mag] <- 1L

  # recurrence pval calculation
  mat_recur_pval <- apply(mat_score_discrete, 2, recur_pval)
  return(mat_recur_pval)
}

#' Get the Recurrence P Values
#'
#' This employs a permutation approach to the recurrence P values. See References for details.
#'
#' @param mat_value_clustered A matrix of expression/alteration clustered by the sparse hierarchical clustering.
#' @param df_feature A data frame that annotates function scores of features. It must contains columns: Feature, FunctionScore. Function score is either 1 (activating) or -1 (inhibiting). If df_feature = NULL, both 1 (activating) and -1 (inhibiting) are considered for each feature.
#' @return A matrix that contains recurrence P values for each element of the input (mat_value_clustered).
#' @references
#'
#' Li X., et al. (2020) Precision combination therapies from recurrent oncogenic co-alterations \doi{https://doi.org/10.1101/2020.06.03.132514}.
#'
#' @examples
#' library(reflect)
#' mat_value <- egfr_data$mat_value
#' wbound <- 2.0
#' mat_value_clustered <- sparse_hclust(mat_value, wbound)$mat_value_clustered
#'
#' df_feature <- egfr_data$df_feature
#' mat_recur_pval <- get_recur_pval(mat_value_clustered, df_feature)
#' @export
get_recur_pval <- function(mat_value_clustered,
                           df_feature = NULL
) {
  if (sum(duplicated(rownames(mat_value_clustered))) > 0 |
      sum(duplicated(colnames(mat_value_clustered))) > 0) {
    stop('samples and features (rownames and colnames of mat_value_clustered) need to be unique')
  }

  if (!is.null(df_feature)) {
    if (length(setdiff(c('Feature', 'FunctionScore'), colnames(df_feature))) > 0) {
      stop('df_feature must have columns: Feature, FunctionScore')
    }

    if (length(setdiff(colnames(mat_value_clustered), df_feature$Feature)) > 0) {
      stop('mat_value_clustered must has the features contained in df_feature')
    }

    vec_fs <- df_feature$FunctionScore[match(colnames(mat_value_clustered), df_feature$Feature)]
    mat_recur_pval <- get_recur_pval_kernel(mat_value_clustered, vec_fs)
  } else {
    if (min(mat_value_clustered) < 0 & max(mat_value_clustered) > 0) {
      vec_fs <- rep(1, ncol(mat_value_clustered))
      mat_recur_pval_gain <- get_recur_pval_kernel(mat_value_clustered, vec_fs)
      vec_fs <- rep(-1, ncol(mat_value_clustered))
      mat_recur_pval_loss <- get_recur_pval_kernel(mat_value_clustered, vec_fs)
      mat_recur_pval <- mat_recur_pval_gain * mat_recur_pval_loss
    } else if (min(mat_value_clustered) < 0) {
      vec_fs <- rep(-1, ncol(mat_value_clustered))
      mat_recur_pval <- get_recur_pval_kernel(mat_value_clustered, vec_fs)
    } else if (max(mat_value_clustered) > 0) {
      vec_fs <- rep(1, ncol(mat_value_clustered))
      mat_recur_pval <- get_recur_pval_kernel(mat_value_clustered, vec_fs)
    }
  }

  return(mat_recur_pval)
}


# get recurrent/actionable features
#' Get the Recurrent and Actionable Features for Each Sample
#'
#' Given a matrix of expression/alteration, a matrix of recurrence P values, a P value threshold, and optionally a set of actionable features, filter the recurrent and actionable features for each sample from the cohort.
#'
#' @param mat_value A matrix of expression/alteration with samples as rows and features as columns.
#' @param mat_recur_pval A matrix that contains recurrence P values. It may not have the same row/column orders as mat_value, but the set of samples/features must be same between the two matrices.
#' @param pval_threshold The threshold of P value below which are considered statistically significant. Default 0.05.
#' @param df_feature A data frame that annotates actionability of features. It must contains columns: Feature, IsActionable. The type of IsActionable is logical. If df_feature = NULL, all features are considered to be actionable.
#' @return \item{df_recur_actionable}{A data frame of samples having recurrent and actionable features. It contains columns: SampleID, Feature, Feature_Value, Feature_Recur_Pval.} \item{mat_value_recur_actionable}{A matrix that only has recurrent and actionable element as in the input (mat_value).}
#' @examples
#' library(reflect)
#' mat_value <- egfr_data$mat_value
#' wbound <- 2.0
#' mat_value_clustered <- sparse_hclust(mat_value, wbound)$mat_value_clustered
#'
#' df_feature <- egfr_data$df_feature
#' mat_recur_pval <- get_recur_pval(mat_value_clustered, df_feature)
#'
#' recur_actionable <- get_recur_actionable_features(mat_value, mat_recur_pval)
#' @export
get_recur_actionable_features <- function(mat_value,
                                              mat_recur_pval,
                                              pval_threshold = 0.05,
                                              df_feature = NULL
) {
  if (sum(duplicated(rownames(mat_value))) > 0 |
      sum(duplicated(colnames(mat_value))) > 0) {
    stop('samples and features (rownames and colnames of mat_value) need to be unique')
  }

  # test if mat_value and mat_recur_pval have same samples
  if (!all.equal(sort(rownames(mat_value)), sort(rownames(mat_recur_pval)))) {
    stop('mat_value and mat_recur_pval must have the same samples')
  }

  # test if features in mat_recur_pval are in mat_value
  if (length(setdiff(colnames(mat_recur_pval), colnames(mat_value))) > 0) {
    stop('features in mat_recur_pval must be in mat_value')
  }

  # test actionable features
  if (!is.null(df_feature)) {
    if (length(setdiff(c('Feature', 'IsActionable'), colnames(df_feature))) > 0) {
      stop('df_feature must have columns: Feature, IsActionable')
    }

    if (length(setdiff(colnames(mat_value), df_feature$Feature)) > 0) {
      stop('mat_value_clustered must has the features contained in df_feature')
    }
  }

  # convert to tibble and combine
  tb_data <- dplyr::tibble(SampleID = rownames(mat_value)) %>%
    dplyr::bind_cols(dplyr::as_tibble(mat_value)) %>%
    tidyr::gather(Feature, Feature_Value, -SampleID)
  tb_recur <- dplyr::tibble(SampleID = rownames(mat_recur_pval)) %>%
    dplyr::bind_cols(dplyr::as_tibble(mat_recur_pval)) %>%
    tidyr::gather(Feature, Feature_Recur_Pval, -SampleID)
  tb <- tb_data %>%
    dplyr::inner_join(tb_recur, by = c('SampleID', 'Feature')) %>%
    dplyr::arrange(SampleID, Feature_Recur_Pval) %>%
    dplyr::select(SampleID, Feature, Feature_Value, Feature_Recur_Pval)

  # filter p values
  tb <- tb %>%
    dplyr::filter(Feature_Recur_Pval <= pval_threshold)

  # filter actionable
  if (!is.null(df_feature)) {
    tb <- tb %>%
      dplyr::filter(Feature %in% df_feature$Feature[df_feature$IsActionable])
  }
  tb_recur_actionable <- tb

  # matrix of score with recurrent and actionable features
  tb <- tb %>%
    dplyr::select(SampleID, Feature, Feature_Value) %>%
    tidyr::spread(key = Feature, value = Feature_Value)
  tb <- dplyr::tibble(SampleID = rownames(mat_recur_pval)) %>%
    dplyr::left_join(tb, by = 'SampleID') %>%
    replace(is.na(.), 0)
  mat_value_recur_actionable <- as.matrix(tb[, -1])
  rownames(mat_value_recur_actionable) <- tb$SampleID
  valid_features <- colnames(mat_recur_pval)[colnames(mat_recur_pval) %in% colnames(mat_value_recur_actionable)]
  mat_value_recur_actionable <- mat_value_recur_actionable[, valid_features]

  return(list('df_recur_actionable' = as.data.frame(tb_recur_actionable),
              'mat_value_recur_actionable' = mat_value_recur_actionable))
}


# get co-altered targets
#' Get Co-altered Targets from Stratification Status and Recurrent/Actionable Features
#'
#' For each sample, given its stratification status and recurrent/actionable features identified by REFLECT, it combines them to form co-altered targets towards combination therapies.
#'
#' @param df_sample A data frame of samples including tumor types and stratification status. It must contain columns: SampleID, TumorType, Stratification.
#' @param df_recur_actionable A data frame of samples having recurrent and actionable features. It must contain columns: SampleID, Feature, Feature_Value, Feature_Recur_Pval. It can be generated by function get_recur_actionable_features.
#' @return A data frame of samples that have both stratification status and recurrent/actionable features for each sample. It contains columns: SampleID, TumorType, Stratification, Feature, Feature_Value, Feature_Recur_Pval.
#' @examples
#' library(reflect)
#' mat_value <- egfr_data$mat_value
#' wbound <- 2.0
#' mat_value_clustered <- sparse_hclust(mat_value, wbound)$mat_value_clustered
#'
#' df_feature <- egfr_data$df_feature
#' mat_recur_pval <- get_recur_pval(mat_value_clustered, df_feature)
#'
#' recur_actionable <- get_recur_actionable_features(mat_value, mat_recur_pval)
#' df_recur_actionable <- recur_actionable$df_recur_actionable
#'
#' df_sample <- egfr_data$df_sample
#' df_coaltered_targets <- get_coaltered_targets(df_sample, df_recur_actionable)
#' @export
get_coaltered_targets <- function(df_sample,
                                  df_recur_actionable
) {
  if (length(setdiff(c('SampleID', 'TumorType', 'Stratification'), colnames(df_sample))) > 0) {
    stop('df_sample must have columns: SampleID, TumorType, Stratification')
  }

  # convert factors to strings
  if (is.factor(df_sample$SampleID)) {
    df_sample$SampleID <- as.character(df_sample$SampleID)
  }
  if (is.factor(df_sample$TumorType)) {
    df_sample$TumorType <- as.character(df_sample$TumorType)
  }
  if (is.factor(df_sample$Stratification)) {
    df_sample$Stratification <- as.character(df_sample$Stratification)
  }

  tb <- dplyr::as_tibble(df_sample) %>%
    dplyr::select(SampleID, TumorType, Stratification)

  # convert factors to strings
  if (is.factor(df_recur_actionable$SampleID)) {
    df_recur_actionable$SampleID <- as.character(df_recur_actionable$SampleID)
  }

  tb <- dplyr::as_tibble(df_recur_actionable) %>%
    dplyr::left_join(tb, by = "SampleID") %>%
    dplyr::select(SampleID, TumorType, Stratification, everything())

  df_coaltered_targets <- as.data.frame(tb)
  return(df_coaltered_targets)
}


# pipelines ------
#' Run an End-to-end REFLECT Pipeline
#'
#' Given an expression/alteration matrix and annotations of samples and features, the REFLECT pipeline nominates co-altered, recurrent, actionable combination targets. The tunning parameter wbound is selected based on a permuation approach. Note this pipeline may take a long time because obtaining a tunning parameter through gap statitic may be computationally expensive.
#'
#' @param mat_value A matrix of expression/alteration with samples as rows and features as columns.
#' @param df_sample A data frame of samples including tumor types and stratification status. It must contain columns: SampleID, TumorType, Stratification.
#' @param df_feature A data frame that annotates function scores and actionabilities of features. It must contains columns: Feature, FunctionScore, IsActionable. If df_feature = NULL, both 1 (activating) and -1 (inhibiting) are considered as function scores for each feature, and all features are considered as being actionable.
#' @param dissimilarity A string for the type of dissimilarity, either "squared.distance" or "absolute.value". Default "squared.distance".
#' @param wbounds The sequence of tuning parameters to consider. If NULL, then a default sequence seq(1.1, sqrt(ncol(mat_value)), 100) will be used. If non-null, should be greater than 1.
#' @param nperms The number of permutations to perform. Default 10.
#' @param min_number_features The minimal number of features that the best wbound could generate. Only wbounds that generates more than this number of features can considered. Default 10.
#' @param pval_threshold The threshold of P value below which are considered statistically significant. Default 0.05.
#' @return \item{gapstat_bestwbound}{Gap statistic profile and the value of best tunning parameter wbound. See detial in function get_best_wbound.} \item{shc}{An object of class hclust which describes the tree produced by the clustering process. See detial in function sparse_hclust.} \item{mat_recur_pval}{A matrix of recurrence P values. See detial in function get_recur_pval.} \item{recur_actionable}{An object for rcurrent and acitionable features. See detial in function get_recur_actionable_features.} \item{df_coaltered_targets}{A data frame of samples that have both stratification status and recurrent/actionable features for each sample. See detial in function get_coaltered_targets.}
#' @examples
#' library(reflect)
#' mat_value <- egfr_data$mat_value
#' df_sample <- egfr_data$df_sample
#' df_feature <- egfr_data$df_feature
#' res <- reflect_pipeline(mat_value, df_sample, df_feature)
#' @export
reflect_pipeline <- function(mat_value,
                             df_sample,
                             df_feature = NULL,
                             dissimilarity = c("squared.distance",
                                               "absolute.value"),
                             wbounds = NULL,
                             nperms = 10,
                             min_number_features = 10,
                             pval_threshold = 0.05) {
  if (sum(duplicated(rownames(mat_value))) > 0 |
      sum(duplicated(colnames(mat_value))) > 0) {
    stop('samples and features (rownames and colnames of mat_value) need to be unique')
  }

  if (length(setdiff(c('SampleID', 'TumorType', 'Stratification'), colnames(df_sample))) > 0) {
    stop('df_sample must have columns: SampleID, TumorType, Stratification')
  }
  if (length(setdiff(rownames(mat_value), df_sample$SampleID)) > 0) {
    stop('mat_value must has the samples contained in df_sample')
  }

  if (!is.null(df_feature)) {
    if (length(setdiff(c('Feature', 'FunctionScore'), colnames(df_feature))) > 0) {
      stop('df_feature must have columns: Feature, FunctionScore')
    }
    if (length(setdiff(colnames(mat_value), df_feature$Feature)) > 0) {
      stop('mat_value must has the features contained in df_feature')
    }
  }

  gapstat_bestwbound <- get_best_wbound(mat_value,
                                   dissimilarity = dissimilarity,
                                   wbounds = wbounds,
                                   nperms = nperms,
                                   min_number_features = min_number_features)
  best_wbound <- gapstat_bestwbound$best_wbound

  shc <- sparse_hclust(mat_value,
                       best_wbound,
                       dissimilarity = dissimilarity)

  mat_recur_pval <- get_recur_pval(shc$mat_value_clustered,
                                   df_feature = df_feature)

  recur_actionable <-
    get_recur_actionable_features(mat_value,
                                      mat_recur_pval,
                                      pval_threshold,
                                      df_feature)

  df_coaltered_targets <- get_coaltered_targets(df_sample,
                                                recur_actionable$df_recur_actionable)

  return(list("gapstat_bestwbound" = gapstat_bestwbound,
              "shc" = shc,
              "mat_recur_pval" = mat_recur_pval,
              "recur_actionable" = recur_actionable,
              "df_coaltered_targets" = df_coaltered_targets))
}

#' Run a REFLECT Pipeline Given a Precomputed Tunning Parameter
#'
#' Given a precomputed tunning parameter, an expression/alteration matrix and annotations of samples and features, the REFLECT pipeline nominates co-altered, recurrent, actionable combination targets. The tunning parameter wbound is selected based on a permuation approach. Since obtaining a tunning parameter through gap statitic may take a long time, this pipeline can same computational cost if the tunning parameter is precomputed.
#'
#' @param wbound The precalculated tunning parameter for sparse hierarchial clustering.
#' @param mat_value A matrix of expression/alteration with samples as rows and features as columns.
#' @param df_sample A data frame of samples including tumor types and stratification status. It must contain columns: SampleID, TumorType, Stratification.
#' @param df_feature A data frame that annotates function scores and actionabilities of features. It must contains columns: Feature, FunctionScore, IsActionable. If df_feature = NULL, both 1 (activating) and -1 (inhibiting) are considered as function scores for each feature, and all features are considered as being actionable.
#' @param dissimilarity A string for the type of dissimilarity, either "squared.distance" or "absolute.value". Default "squared.distance".
#' @param nperms The number of permutations to perform. Default 10.
#' @param min_number_features The minimal number of features that the best wbound could generate. Only wbounds that generates more than this number of features can considered. Default 10.
#' @param pval_threshold The threshold of P value below which are considered statistically significant. Default 0.05.
#' @return \item{shc}{An object of class hclust which describes the tree produced by the clustering process. See detial in function sparse_hclust.} \item{mat_recur_pval}{A matrix of recurrence P values. See detial in function get_recur_pval.} \item{recur_actionable}{An object for rcurrent and acitionable features. See detial in function get_recur_actionable_features.} \item{df_coaltered_targets}{A data frame of samples that have both stratification status and recurrent/actionable features for each sample. See detial in function get_coaltered_targets.}
#' @examples
#' library(reflect)
#' mat_value <- egfr_data$mat_value
#' gapstat_bestwbound <- get_best_wbound(mat_value)
#' wbound <- gapstat_bestwbound$best_wbound
#'
#' df_sample <- egfr_data$df_sample
#' df_feature <- egfr_data$df_feature
#' res <- reflect_pipeline2(wbound, mat_value, df_sample, df_feature)
#' @export
reflect_pipeline2 <- function(wbound,
                              mat_value,
                              df_sample,
                              df_feature = NULL,
                              dissimilarity = c("squared.distance",
                                                "absolute.value"),
                              nperms = 10,
                              min_number_features = 10,
                              pval_threshold = 0.05) {
  if (sum(duplicated(rownames(mat_value))) > 0 |
      sum(duplicated(colnames(mat_value))) > 0) {
    stop('samples and features (rownames and colnames of mat_value) need to be unique')
  }

  if (length(setdiff(c('SampleID', 'TumorType', 'Stratification'), colnames(df_sample))) > 0) {
    stop('df_sample must have columns: SampleID, TumorType, Stratification')
  }
  if (length(setdiff(rownames(mat_value), df_sample$SampleID)) > 0) {
    stop('mat_value must has the samples contained in df_sample')
  }

  if (!is.null(df_feature)) {
    if (length(setdiff(c('Feature', 'FunctionScore'), colnames(df_feature))) > 0) {
      stop('df_feature must have columns: Feature, FunctionScore')
    }
    if (length(setdiff(colnames(mat_value), df_feature$Feature)) > 0) {
      stop('mat_value must has the features contained in df_feature')
    }
  }

  shc <- sparse_hclust(mat_value,
                       wbound,
                       dissimilarity = dissimilarity)
  mat_recur_pval <- get_recur_pval(shc$mat_value_clustered,
                                   df_feature = df_feature)
  recur_actionable <-
    get_recur_actionable_features(mat_value,
                                      mat_recur_pval,
                                      pval_threshold,
                                      df_feature)

  df_coaltered_targets <- get_coaltered_targets(df_sample,
                                                recur_actionable$df_recur_actionable)

  return(list("shc" = shc,
              "mat_recur_pval" = mat_recur_pval,
              "recur_actionable" = recur_actionable,
              "df_coaltered_targets" = df_coaltered_targets))
}
