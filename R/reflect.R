# bestw function -------
#' @export
get_gapstat_bestw <- function(mat, 
                              dissimilarity = c("squared.distance",
                                                "absolute.value"), 
                              nwbounds = 100, 
                              nperms = 10,
                              min_number_features = 10
) {
  if (sum(duplicated(rownames(mat))) > 0 |
      sum(duplicated(colnames(mat))) > 0) {
    stop('samples and features (rownames and colnames of mat) need to be unique')
  }
  
  dissimilarity <- match.arg(dissimilarity)
  
  wbounds <- seq(1.1, sqrt(ncol(mat)) * 1.0, len = nwbounds)
  perm.out <- sparcl::HierarchicalSparseCluster.permute(mat, 
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
  
  # get bestw ------
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
  bestw <- tb$Wbound[1]
  
  return(list("df_gapstat" = df_gapstat,
              "bestw" = bestw))
}


# sparse_hclust -----
#' @export
sparse_hclust <- function(mat,
                          bestw,
                          dissimilarity = c("squared.distance",
                                            "absolute.value")
) {
  if (sum(duplicated(rownames(mat))) > 0 |
      sum(duplicated(colnames(mat))) > 0) {
    stop('samples and features (rownames and colnames of mat) need to be unique')
  }
  
  sparsehc <- sparcl::HierarchicalSparseCluster(x = mat,
                                                dissimilarity = dissimilarity,
                                                dists = NULL,
                                                wbound = bestw)
  
  # the original hc has numbers as lables
  # change the labels to sampleid
  sparsehc$hc$labels <- rownames(mat)
  
  weight <- sparsehc$ws
  names(weight) <- colnames(mat)
  weight <- sort(weight, decreasing = TRUE)
  
  features_nonzero <- names(weight)[weight > 0]
  mat_clustered <- mat[sparsehc$hc$order, features_nonzero]
  
  return(list("hc" = sparsehc$hc, 
              "weight" = weight, 
              "features_nonzero" = features_nonzero,
              "mat_clustered" = mat_clustered))
}


# robustness test -------
get_top_n_features <- function(mat, n_top,
                               dissimilarity, nwbounds, nperms,
                               min_number_features
) {
  gapstat_bestw <- get_gapstat_bestw(mat,
                                     dissimilarity = dissimilarity, 
                                     nwbounds = nwbounds, 
                                     nperms = nperms,
                                     min_number_features = min_number_features) 
  bestw <- gapstat_bestw$bestw
  sparsehc <- sparcl::HierarchicalSparseCluster(x = mat,
                                                dissimilarity = dissimilarity,
                                                dists = NULL,
                                                wbound = bestw)
  tb <- dplyr::tibble(Feature = colnames(mat), ws = sparsehc$ws[, 1])
  tb <- tb %>%
    dplyr::filter(ws > 0) %>%
    dplyr::arrange(dplyr::desc(ws))
  if(nrow(tb) > n_top) tb <- tb[1:n_top, ]
  return(tb$Feature)
}

get_score <- function(mat, indices, features_ref, n_top_sub,
                      dissimilarity, nwbounds, nperms,
                      min_number_features){
  mat_sub <- mat[indices, ]
  features_sub <- get_top_n_features(mat_sub, n_top_sub,
                                     dissimilarity, nwbounds, nperms,
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

#' @export
subsample_robustness <- function(mat,
                                 dissimilarity = c("squared.distance",
                                                   "absolute.value"), 
                                 nwbounds = 100, 
                                 nperms = 10,
                                 min_number_features = 10,
                                 n_top_ref = 20, 
                                 n_top_sub = 10,
                                 fs_subsample = seq(0.6, 0.95, 0.05),
                                 n_iters = 200,
                                 no_cores = NULL,
                                 seed = 123
) {
  if (sum(duplicated(rownames(mat))) > 0 |
      sum(duplicated(colnames(mat))) > 0) {
    stop('samples and features (rownames and colnames of mat) need to be unique')
  }
  
  dissimilarity <- match.arg(dissimilarity)
  
  # create reference for n_top_ref from the full dataset
  features_ref <- get_top_n_features(mat, n_top_ref,
                                     dissimilarity, nwbounds, nperms,
                                     min_number_features)
  
  # subsampling and calculate robustness score ------
  n_totsample <- nrow(mat)
  
  indices_list <- list()
  f_subsample_vec <- c()
  set.seed(seed)
  for (f_subsample in fs_subsample) {
    for (i in c(1:n_iters)) {
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
                          varlist=c("get_gapstat_bestw", "get_top_n_features", "get_score",
                                    "indices_list", "mat", "features_ref", "n_top_sub",
                                    "dissimilarity", "nwbounds", "nperms",
                                    "min_number_features"), 
                          envir=environment())
  
  get_score_wrap <- function(id){
    indices <- indices_list[[id]]
    get_score(mat, indices, features_ref, n_top_sub, 
              dissimilarity, nwbounds, nperms,
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
get_recur_pval_kernel <- function(mat, vec_fs) {
  # score matrix
  mat_fs <- matrix(rep(vec_fs, nrow(mat)),
                   ncol = length(vec_fs), byrow = TRUE)
  colnames(mat_fs) <- names(mat)
  rownames(mat_fs) <- rownames(mat)
  mat_score <- mat * mat_fs
  
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

#' @export
get_recur_pval <- function(mat, 
                           df_feature = NULL
) {
  if (sum(duplicated(rownames(mat))) > 0 |
      sum(duplicated(colnames(mat))) > 0) {
    stop('samples and features (rownames and colnames of mat) need to be unique')
  }
  
  if (!is.null(df_feature)) {
    if (length(setdiff(c('Feature', 'FunctionScore'), colnames(df_feature))) > 0) {
      stop('df_feature must have columns: Feature, FunctionScore')
    }
    
    if (length(setdiff(colnames(mat), df_feature$Feature)) > 0) {
      stop('mat must has the features contained in df_feature')
    }
    
    vec_fs <- df_feature$FunctionScore[match(colnames(mat), df_feature$Feature)]
    mat_recur_pval <- get_recur_pval_kernel(mat, vec_fs)
  } else {
    if (min(mat) < 0 & max(mat) > 0) {
      vec_fs <- rep(1, ncol(mat))
      mat_recur_pval_gain <- get_recur_pval_kernel(mat, vec_fs)
      vec_fs <- rep(-1, ncol(mat))
      mat_recur_pval_loss <- get_recur_pval_kernel(mat, vec_fs)
      mat_recur_pval <- mat_recur_pval_gain * mat_recur_pval_loss
    } else if (min(mat) < 0) {
      vec_fs <- rep(-1, ncol(mat))
      mat_recur_pval <- get_recur_pval_kernel(mat, vec_fs)
    } else if (max(mat) > 0) {
      vec_fs <- rep(1, ncol(mat))
      mat_recur_pval <- get_recur_pval_kernel(mat, vec_fs)
    }
  }
  
  return(mat_recur_pval)
}


# get recurrent/actionable features
#' @export
get_recurrent_actionable_features <- function(mat, 
                                              mat_recur_pval,
                                              pval_threshold = 0.05,
                                              actionable_features = NULL
) {
  if (sum(duplicated(rownames(mat))) > 0 |
      sum(duplicated(colnames(mat))) > 0) {
    stop('samples and features (rownames and colnames of mat) need to be unique')
  }
  
  # test if mat and mat_recur_pval have same samples
  if (!all.equal(sort(rownames(mat)), sort(rownames(mat_recur_pval)))) {
    stop('mat and mat_recur_pval must have the same samples')
  }
  
  # test if features in mat_recur_pval are in mat
  if (length(setdiff(colnames(mat_recur_pval), colnames(mat))) > 0) {
    stop('features in mat_recur_pval must be in mat')
  }
  
  # test actionable features
  if (!is.null(actionable_features)) {
    if (length(setdiff(actionable_features, colnames(mat))) > 0) {
      stop('some features in actionable_features are not contained in mat')
    }
  }
  
  # convert to tibble and combine
  tb_data <- dplyr::tibble(SampleID = rownames(mat)) %>%
    dplyr::bind_cols(dplyr::as_tibble(mat)) %>%
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
  if (!is.null(actionable_features)) {
    tb <- tb %>%
      dplyr::filter(Feature %in% actionable_features)
  }
  tb_recur_actionable <- tb
  
  # matrix of score with recurrent and actionable features
  tb <- tb %>%
    dplyr::select(SampleID, Feature, Feature_Value) %>%
    tidyr::spread(key = Feature, value = Feature_Value)
  tb <- dplyr::tibble(SampleID = rownames(mat_recur_pval)) %>%
    dplyr::left_join(tb, by = 'SampleID') %>%
    replace(is.na(.), 0)
  mat_recur_actionable <- as.matrix(tb[, -1])
  rownames(mat_recur_actionable) <- tb$SampleID
  valid_features <- colnames(mat_recur_pval)[colnames(mat_recur_pval) %in% colnames(mat_recur_actionable)]
  mat_recur_actionable <- mat_recur_actionable[, valid_features]

  return(list('df_recur_actionable' = as.data.frame(tb_recur_actionable),
              'mat_recur_actionable' = mat_recur_actionable))
}


# get co-altered targets
#' @export
get_coaltered_targets <- function(df_sample,
                                  df_recur_actionable
) {
  if (length(setdiff(c('SampleID', 'TumorType', 'Stratification'), colnames(df_sample))) > 0) {
    stop('df_sample must have columns: SampleID, TumorType, Stratification')
  }
  
  tb <- dplyr::as_tibble(df_sample) %>%
    dplyr::select(SampleID, TumorType, Stratification)
  
  tb <- dplyr::as_tibble(df_recur_actionable) %>%
    dplyr::left_join(tb, by = "SampleID") %>%
    dplyr::select(SampleID, TumorType, Stratification, everything())
  
  df_coaltered_targets <- as.data.frame(tb)
  return(df_coaltered_targets)
}


# pipelines ------
#' @export
reflect_pipeline <- function(mat,
                             df_sample,
                             df_feature = NULL,
                             actionable_features = NULL,
                             dissimilarity = c("squared.distance",
                                               "absolute.value"),
                             nwbounds = 100, 
                             nperms = 10,
                             min_number_features = 10,
                             pval_threshold = 0.05) {
  if (sum(duplicated(rownames(mat))) > 0 |
      sum(duplicated(colnames(mat))) > 0) {
    stop('samples and features (rownames and colnames of mat) need to be unique')
  }
  
  if (length(setdiff(c('SampleID', 'TumorType', 'Stratification'), colnames(df_sample))) > 0) {
    stop('df_sample must have columns: SampleID, TumorType, Stratification')
  }
  if (length(setdiff(rownames(mat), df_sample$SampleID)) > 0) {
    stop('mat must has the samples contained in df_sample')
  }
  
  if (!is.null(df_feature)) {
    if (length(setdiff(c('Feature', 'FunctionScore'), colnames(df_feature))) > 0) {
      stop('df_feature must have columns: Feature, FunctionScore')
    }
    if (length(setdiff(colnames(mat), df_feature$Feature)) > 0) {
      stop('mat must has the features contained in df_feature')
    }
  }
  
  gapstat_bestw <- get_gapstat_bestw(mat, 
                                     dissimilarity = dissimilarity, 
                                     nwbounds = nwbounds, 
                                     nperms = nperms,
                                     min_number_features = min_number_features)
  bestw <- gapstat_bestw$bestw
  
  shc <- sparse_hclust(mat, 
                       bestw,
                       dissimilarity = dissimilarity)
  
  mat_recur_pval <- get_recur_pval(shc$mat_clustered,
                                   df_feature = df_feature)
  
  recur_actionable <- 
    get_recurrent_actionable_features(mat,
                                      mat_recur_pval,
                                      pval_threshold,
                                      actionable_features)
  
  df_coaltered_targets <- get_coaltered_targets(df_sample,
                                                recur_actionable$df_recur_actionable)
  
  return(list("gapstat_bestw" = gapstat_bestw,
              "shc" = shc,
              "mat_recur_pval" = mat_recur_pval,
              "recur_actionable" = recur_actionable,
              "df_coaltered_targets" = df_coaltered_targets))
}

#' @export
reflect_pipeline2 <- function(bestw,
                              mat,
                              df_sample,
                              df_feature = NULL,
                              actionable_features = NULL,
                              dissimilarity = c("squared.distance",
                                                "absolute.value"),
                              nwbounds = 100, 
                              nperms = 10,
                              min_number_features = 10,
                              pval_threshold = 0.05) {
  if (sum(duplicated(rownames(mat))) > 0 |
      sum(duplicated(colnames(mat))) > 0) {
    stop('samples and features (rownames and colnames of mat) need to be unique')
  }
  
  if (length(setdiff(c('SampleID', 'TumorType', 'Stratification'), colnames(df_sample))) > 0) {
    stop('df_sample must have columns: SampleID, TumorType, Stratification')
  }
  if (length(setdiff(rownames(mat), df_sample$SampleID)) > 0) {
    stop('mat must has the samples contained in df_sample')
  }
  
  if (!is.null(df_feature)) {
    if (length(setdiff(c('Feature', 'FunctionScore'), colnames(df_feature))) > 0) {
      stop('df_feature must have columns: Feature, FunctionScore')
    }
    if (length(setdiff(colnames(mat), df_feature$Feature)) > 0) {
      stop('mat must has the features contained in df_feature')
    }
  }
  
  shc <- sparse_hclust(mat, 
                       bestw,
                       dissimilarity = dissimilarity)
  mat_recur_pval <- get_recur_pval(shc$mat_clustered,
                                   df_feature = df_feature)
  recur_actionable <- 
    get_recurrent_actionable_features(mat,
                                      mat_recur_pval,
                                      pval_threshold,
                                      actionable_features)
  
  df_coaltered_targets <- get_coaltered_targets(df_sample,
                                                recur_actionable$df_recur_actionable)
  
  return(list("shc" = shc,
              "mat_recur_pval" = mat_recur_pval,
              "recur_actionable" = recur_actionable,
              "df_coaltered_targets" = df_coaltered_targets))
}
