# bar chart of sample tumor types --------
#' Plot a Bar Chart of Tumor Type Distribution
#'
#' This function takes a data frame of samples and corresponding tumor types, calculates occurrences of each tumor type, then plots the occurences as a bar chart.
#'
#' @param df A data frame that has columns: SampleID, TumorType.
#' @param main_label Main lable. Default NULL.
#' @examples
#' library(reflect)
#' df <- egfr_data$df_sample[, c("SampleID", "TumorType")]
#' plot_bar_tumortypes(df)
#' @export
plot_bar_tumortypes <- function(df, main_label = NULL) {
  if (nrow(df) == 0) return(NULL)
  if (length(setdiff(c('SampleID', 'TumorType'), colnames(df))) > 0) {
    stop('df must have columns: SampleID, TumorType')
  }

  # convert factors to strings
  if (is.factor(df$SampleID)) {
    df$SampleID <- as.character(df$SampleID)
  }
  if (is.factor(df$TumorType)) {
    df$TumorType <- as.character(df$TumorType)
  }

  tb <- dplyr::as_tibble(df)
  tb <- tb %>%
    dplyr::group_by(TumorType) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::arrange(dplyr::desc(Count))

  # collapse small fractions
  tot_n <- sum(tb$Count)
  tb <- tb %>%
    dplyr::mutate(Fract = 1.0 * Count / tot_n)

  # plot
  main = ifelse(is.null(main_label),
                paste0(sum(tb$Count), " Samples"),
                paste0(main_label, " (", sum(tb$Count), " Samples)"))
  p <- ggpubr::ggbarplot(tb, "TumorType", "Count",
                         fill = "steelblue", color = "steelblue",
                         # label = TRUE, lab.pos = "in", lab.col = "white",
                         title = main) +
    ggpubr::rotate_x_text(angle = 90)
  p
}


#' Plot a Stacked Bar Chart of Tumor Types for TCGA and Cell Line Samples
#'
#' This function takes a data frame of samples and corresponding tumor types, separate TCGA and cell line samples, calculates occurrences of each tumor type, then plots the occurences as a stacked bar chart.
#'
#' @param df A data frame that has columns: SampleID, TumorType.
#' @examples
#' library(reflect)
#' df <- egfr_data$df_sample[, c("SampleID", "TumorType")]
#' plot_bar_tumortypes_stack_tcga_ccl(df)
#' @export
plot_bar_tumortypes_stack_tcga_ccl <- function(df, main_label = NULL) {
  if (nrow(df) == 0) return(NULL)
  if (length(setdiff(c('SampleID', 'TumorType'), colnames(df))) > 0) {
    stop('df must have columns: SampleID, TumorType')
  }

  # convert factors to strings
  if (is.factor(df$SampleID)) {
    df$SampleID <- as.character(df$SampleID)
  }
  if (is.factor(df$TumorType)) {
    df$TumorType <- as.character(df$TumorType)
  }

  tb <- dplyr::as_tibble(df) %>%
    dplyr::select(SampleID, TumorType)
  tb_tcga <- tb %>%
    dplyr::filter(grepl('^TCGA', SampleID))
  tb_ccl <- tb %>%
    dplyr::filter(!grepl('^TCGA', SampleID))

  tb <- dplyr::tibble()
  # tcga
  if (nrow(tb_tcga) > 0) {
    tb_tcga <- tb_tcga %>%
      dplyr::group_by(TumorType) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(Count))
    tb <- tb %>%
      dplyr::bind_rows(tb_tcga)
  }

  # ccl
  if (nrow(tb_ccl) > 0) {
    tb_ccl <- tb_ccl %>%
      dplyr::group_by(TumorType) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(Count))
    tb <- tb %>%
      dplyr::bind_rows(tb_ccl)
  }

  # collaps small fractions
  tot_n <- sum(tb$Count)
  tb <- tb %>%
    dplyr::mutate(Fract = 1.0 * Count / tot_n)

  # plot
  main = ifelse(is.null(main_label),
                paste0(sum(tb$Count), " Samples"),
                paste0(main_label, " (", sum(tb$Count), " Samples)"))
  p <- ggpubr::ggbarplot(tb, "TumorType", "Count",
                         fill = "steelblue", color = "steelblue",
                         # label = TRUE, lab.pos = "in", lab.col = "white",
                         title = main) +
    ggpubr::rotate_x_text(angle = 90)
  p
}


# gap statistic --------
#' Plot a Gap Statistic Profile
#'
#' This function plots a gap statistic profile as a function of the number of non-zero features, and optionally as a function of tunning parameters.
#'
#' @param df A data frame that has columns: Wbound, NumberFeatures, GapStat, SD.
#' @param plot_tuning_parameter A logical value. If TRUE, it also draws a second plot of gap statistic as a function of tunning parameters. Default FALSE.
#' @examples
#' library(reflect)
#' df <- egfr_result$gapstat_bestwbound$df_gapstat
#' plot_gapstat(df)
#' @export
plot_gapstat <- function(df, plot_tuning_parameter = FALSE) {
  if (length(setdiff(c("Wbound", "NumberFeatures", "GapStat", "SD"), colnames(df))) > 0) {
    stop('df must have columns: Wbound, NumberFeatures, GapStat, SD')
  }

  tb <- dplyr::as_tibble(df)

  # keep maximum gap statistic for each number of features
  tb <- tb %>%
    dplyr::arrange(NumberFeatures, dplyr::desc(GapStat)) %>%
    dplyr::distinct(NumberFeatures, .keep_all = TRUE)

  p1 <- ggpubr::ggline(tb, x = "NumberFeatures", y = "GapStat",
                       numeric.x.axis = TRUE,
                       xlab = "Number of Features",
                       ylab = "Gap Statistic") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = GapStat - SD, ymax = GapStat + SD))

  if (plot_tuning_parameter) {
    p2 <- ggpubr::ggline(tb, x = "Wbound", y = "GapStat",
                         numeric.x.axis = TRUE,
                         xlab = "Tuning Parameter",
                         ylab = "Gap Statistic") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = GapStat - SD, ymax = GapStat + SD))
    ggpubr::ggarrange(p1, p2,
                      ncol = 1, nrow = 2)
  } else {
    p1
  }
}


# bar chart of weights --------
#' Plot a Bar Chart of Feature Weights
#'
#' This function plots weights of non-zero features.
#'
#' @param df A data frame that has columns: Feature, Weight.
#' @param n_top An integer number showing how many top features are displayed. If NULL, all non-zero features are displayed.
#' @examples
#' library(reflect)
#' df <- data.frame(Feature = names(egfr_result$shc$weight),
#'                  Weight = egfr_result$shc$weight)
#' plot_bar_weights(df)
#' @export
plot_bar_weights <- function(df, n_top = NULL) {
  if (length(setdiff(c('Feature', 'Weight'), colnames(df))) > 0) {
    stop('df must have columns: Feature, Weight')
  }

  # convert factors to strings
  if (is.factor(df$Feature)) {
    df$Feature <- as.character(df$Feature)
  }

  tb <- dplyr::as_tibble(df)
  tb <- tb %>%
    dplyr::filter(Weight > 0) %>%
    dplyr::arrange(dplyr::desc(Weight))

  if (!is.null(n_top)) {
    tb <- tb[1:n_top, ]
  }

  # plot
  p <- ggpubr::ggbarplot(tb, "Feature", "Weight",
                         fill = "steelblue", color = "steelblue") +
    ggpubr::rotate_x_text(angle = 90)
  p
}


# heatmap ------
#' Plot a Heatmap with Feature and Sample Covariate Bars
#'
#' This function plots a heatmap of expression/alteration values. In addition, samples' tumor types and categories are plotted as row covariate bars, and features' weights are plotted as a column covariate bar.
#'
#' @param mat_value A matrix of expression/alteration values with samples as rows and features as columns.
#' @param name A string that is the name of the legend. Default 'Value'.
#' @param discrete A logical value. If mat_value is discrete, set it as TRUE. Default FALSE.
#' @param df_sample A data frame of sample annotations. It contains columns: SampleID, TumorType (optional), Category (optional). If NULL, no row covariate bar is plotted.
#' @param df_weight A data frame of feature weights. It must contain columns: Feature, Weight. If NULL, no column covariate bar is plotted.
#' @param order_features A logical value indicating whether features are ordered according to weights. Default TRUE.
#' @param trim_zero_features A logical value indicating whether the zero features are removed. Default TRUE.
#' @param n_top_features An integer showing how many top features are displayed. If NULL, all features are displayed.
#' @param df_tcga_color A data frame of TCGA tumor type colors. It must contain columns: TumorType, Color. If NULL, colors are automatically generated.
#' @param df_ccl_color A data frame of cell line tumor type colors. It must contain columns: TumorType, Color. If NULL, colors are automatically generated.
#' @param df_category_color A data frame of category colors. It must contain columns: Category, Color. If NULL, colors are automatically generated.
#' @param category_name A string that is displayed as the name of the category covariate bar. If NULL, the displayed name is "Category".
#' @param cluster_rows A dendrogram or a logical value. If a dendrogram is provided, rows will be clustered according to it. If FALSE, rows are not to be clustered. If TRUE, rows will be clustered by standard hierarchical clustering.
#' @param cluster_columns A dendrogram or a logical value. If a dendrogram is provided, columns will be clustered according to it. If FALSE, columns are ordered with descreasing weigths. If TRUE, rows will be clustered by standard hierarchical clustering.
#' @param show_row_names A logical value indicating whether row names are displayed. Default TRUE.
#' @param show_column_names A logical value indicating whether column names are displayed. Default TRUE.
#' @examples
#' library(reflect)
#' mat_value <- egfr_data$mat_value
#' df_sample <- egfr_data$df_sample
#' df_weight <- data.frame(Feature = names(egfr_result$shc$weight),
#'                         Weight = egfr_result$shc$weight)
#' cluster_rows <- as.dendrogram(egfr_result$shc$hc)
#' df_tcga_color <- reflect_color$df_tcga_color
#' df_ccl_color <- reflect_color$df_ccl_color
#' df_category_color <- reflect_color$df_category_color
#' plot_heatmap(mat_value,
#'              df_sample = df_sample,
#'              df_weight = df_weight,
#'              cluster_rows = cluster_rows,
#'              df_tcga_color = df_tcga_color,
#'              df_ccl_color = df_ccl_color,
#'              df_category_color = df_category_color,
#'              category_name = "Variants",
#'              show_row_names = FALSE)
#' @export
plot_heatmap <- function(mat_value,
                         name = 'Value',
                         discrete = FALSE,
                         df_sample = NULL,
                         df_weight = NULL,
                         order_features = TRUE,
                         trim_zero_features = TRUE,
                         n_top_features = NULL,
                         df_tcga_color = NULL,
                         df_ccl_color = NULL,
                         df_category_color = NULL,
                         category_name = NULL,
                         cluster_rows = FALSE,
                         cluster_columns = FALSE,
                         show_row_names = TRUE,
                         show_column_names = TRUE
) {
  # test consistency between mat_value and df_sample
  if (!is.null(df_sample)) {
    if (!('SampleID' %in% colnames(df_sample))) {
      stop('df_sample must have column: SampleID')
    } else {
      if (length(setdiff(rownames(mat_value), df_sample$SampleID)) > 0) {
        stop('mat_value must has the samples contained in df_sample')
      }
    }

    # convert factors to strings
    if (is.factor(df_sample$SampleID)) {
      df_sample$SampleID <- as.character(df_sample$SampleID)
    }

    # reorder df_sample
    df_sample <- df_sample[match(rownames(mat_value), df_sample$SampleID), ]
  }

  # flag
  has_tumortype <- (!is.null(df_sample)) & ('TumorType' %in% colnames(df_sample))
  if (has_tumortype) {
    # convert factors to strings
    if (is.factor(df_sample$TumorType)) {
      df_sample$TumorType <- as.character(df_sample$TumorType)
    }
  }
  has_category <- (!is.null(df_sample)) & ('Category' %in% colnames(df_sample))
  if (has_category) {
    # convert factors to strings
    if (is.factor(df_sample$Category)) {
      df_sample$Category <- as.character(df_sample$Category)
    }
  }

  # test consistency between mat_value and df_weight
  if (!is.null(df_weight)) {
    if (!('Feature' %in% colnames(df_weight))) {
      stop('df_weight must have column: Feature')
    } else {
      if (length(setdiff(colnames(mat_value), df_weight$Feature)) > 0) {
        stop('mat_value must has the features contained in df_weight')
      }
    }

    # convert factors to strings
    if (is.factor(df_weight$Feature)) {
      df_weight$Feature <- as.character(df_weight$Feature)
    }

    # reorder df_weight
    df_weight <- df_weight[match(colnames(mat_value), df_weight$Feature), ]
  }

  # flag
  has_weight <- (!is.null(df_weight)) & ('Weight' %in% colnames(df_weight))

  # order or trim features
  if (has_weight) {
    if (order_features) {
      index <- order(df_weight$Weight, decreasing = TRUE)
      mat_value <- mat_value[, index]
      df_weight <- df_weight[index, ]
    }
    if (trim_zero_features) {
      index <- c(1:nrow(df_weight))[df_weight$Weight > 0]
      mat_value <- mat_value[, index]
      df_weight <- df_weight[index, ]
    }
  }

  # select top feature if asked for
  if (!is.null(n_top_features)) {
    mat_value <- mat_value[, c(1:min(n_top_features, ncol(mat_value)))]
    df_weight <- df_weight[c(1:min(n_top_features, nrow(df_weight))), ]
  }

  # color map
  if (discrete) {
    col_map <- circlize::colorRamp2(c(-1, 0, 1),
                                    c("blue", "white", "red"))
  } else {
    col_map <- circlize::colorRamp2(c(min(quantile(mat_value, 0.01), -1), 0, max(quantile(mat_value, 0.99), 1)),
                                    c("blue", "white", "red"))
  }

  # all qualitative palettes from RColorBrewer package.
  # Qualitative palettes are supposed to provide X most distinctive colours each.
  # Of course, mixing them joins into one palette also similar colours,
  # but that's the best I can get (74 colors).
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  # show_row_names flag
  if (has_category) {
    show_row_names_main <- FALSE
    show_row_names_tumortype <- FALSE
    if (show_row_names) {
      show_row_names_cat <- TRUE
    } else {
      show_row_names_cat <- FALSE
    }
  } else if (has_tumortype) {
    show_row_names_main <- FALSE
    if (show_row_names) {
      show_row_names_tumortype <- TRUE
    } else {
      show_row_names_tumortype <- FALSE
    }
  } else {
    if (show_row_names) {
      show_row_names_main <- TRUE
    } else {
      show_row_names_main <- FALSE
    }
  }

  # main heatmap
  if (has_weight) {
    col_map_weight = circlize::colorRamp2(c(0, max(df_weight$Weight)), c("white", "green"))
    column_ha <- ComplexHeatmap::HeatmapAnnotation(Weight = df_weight$Weight, col = list(Weight = col_map_weight))
  } else {
    column_ha <- NULL
  }

  ht <- ComplexHeatmap::Heatmap(mat_value,
                                name = name,
                                top_annotation = column_ha,
                                col = col_map,
                                cluster_rows = cluster_rows,
                                cluster_columns = cluster_columns,
                                #row_dend_width = ggplot2::unit(10, "mm"),
                                show_row_names = show_row_names_main,
                                show_column_names = show_column_names,
                                row_title = "Samples",
                                row_title_side = "left",
                                column_title = "Features",
                                column_title_side = "bottom")

  # disease tumortype side heatmap
  if (has_tumortype) {
    tb_tumortype <- dplyr::as_tibble(df_sample) %>%
      dplyr::select(SampleID, TumorType)
    is_tcga <- grepl("^TCGA", tb_tumortype$SampleID)

    # tcga
    if (sum(is_tcga) > 0) {
      tcga_tumortypes <- tb_tumortype %>%
        dplyr::filter(is_tcga) %>%
        .$TumorType %>%
        unique()
      if (is.null(df_tcga_color)) {
        col_map_tcga <- structure(col_vector[1:length(tcga_tumortypes)],
                                  names = tcga_tumortypes)
      } else {
        if (length(setdiff(c('TumorType', 'Color'), colnames(df_tcga_color))) > 0) {
          stop('df_tcga_color must have columns: TumorType, Color')
        }

        # convert factors to strings
        if (is.factor(df_tcga_color$TumorType)) {
          df_tcga_color$TumorType <- as.character(df_tcga_color$TumorType)
        }
        if (is.factor(df_tcga_color$Color)) {
          df_tcga_color$Color <- as.character(df_tcga_color$Color)
        }

        if (length(setdiff(tcga_tumortypes, df_tcga_color$TumorType)) > 0) {
          stop('some TCGA Tumor Types are not included in df_tcga_color')
        }
        col_map_tcga <- structure(df_tcga_color$Color,
                                  names = df_tcga_color$TumorType)
      }
    }

    # ccl color map
    if (sum(!is_tcga) > 0) {
      ccl_tumortypes <- tb_tumortype %>%
        dplyr::filter(!is_tcga) %>%
        .$TumorType %>%
        unique()
      if (is.null(df_ccl_color)) {
        col_map_ccl <- structure(col_vector[1:length(ccl_tumortypes)],
                                 names = ccl_tumortypes)
      } else {
        if (length(setdiff(c('TumorType', 'Color'), colnames(df_ccl_color))) > 0) {
          stop('df_ccl_color must have columns: TumorType, Color')
        }

        # convert factors to strings
        if (is.factor(df_ccl_color$TumorType)) {
          df_ccl_color$TumorType <- as.character(df_ccl_color$TumorType)
        }
        if (is.factor(df_ccl_color$Color)) {
          df_ccl_color$Color <- as.character(df_ccl_color$Color)
        }

        if (length(setdiff(ccl_tumortypes, df_ccl_color$TumorType)) > 0) {
          stop('some CCL Tumor Types are not included in df_ccl_color')
        }
        col_map_ccl <- structure(df_ccl_color$Color,
                                 names = df_ccl_color$TumorType)
      }
    }

    # disease tumortype side heatmap
    if (sum(!is_tcga) == 0) {
      mat_tumortype <- tb_tumortype %>%
        dplyr::select(TumorType) %>%
        as.matrix()
      rownames(mat_tumortype) <- tb_tumortype$SampleID
      mat_tumortype[!is_tcga, 1] <- NA
      colnames(mat_tumortype) <- "TCGA"
      ht_tcga <- ComplexHeatmap::Heatmap(mat_tumortype, name = "TCGA",
                                         col = col_map_tcga,
                                         na_col = "white",
                                         show_row_names = show_row_names_tumortype,
                                         show_column_names = TRUE,
                                         column_title = NULL,
                                         column_title_rot = 90,
                                         column_title_side = "bottom")
      ht_list <- ht + ht_tcga
    } else if (sum(is_tcga) == 0) {
      mat_tumortype <- tb_tumortype %>%
        dplyr::select(TumorType) %>%
        dplyr::as.matrix()
      rownames(mat_tumortype) <- tb_tumortype$SampleID
      mat_tumortype[is_tcga, 1] <- NA
      colnames(mat_tumortype) <- "CCL"
      ht_ccl <- ComplexHeatmap::Heatmap(mat_tumortype, name = "CCL",
                                        col = col_map_ccl,
                                        na_col = "white",
                                        show_row_names = show_row_names_tumortype,
                                        show_column_names = TRUE,
                                        column_title = NULL,
                                        column_title_rot = 90,
                                        column_title_side = "bottom")
      ht_list <- ht + ht_ccl
    } else {
      mat_tumortype <- tb_tumortype %>%
        dplyr::select(TumorType) %>%
        as.matrix()
      rownames(mat_tumortype) <- tb_tumortype$SampleID
      mat_tumortype[!is_tcga, 1] <- NA
      colnames(mat_tumortype) <- "TCGA"
      ht_tcga <- ComplexHeatmap::Heatmap(mat_tumortype, name = "TCGA",
                                         col = col_map_tcga,
                                         na_col = "white",
                                         show_row_names = FALSE,
                                         show_column_names = TRUE,
                                         column_title = NULL,
                                         column_title_rot = 90,
                                         column_title_side = "bottom")

      mat_tumortype <- tb_tumortype %>%
        dplyr::select(TumorType) %>%
        as.matrix()
      rownames(mat_tumortype) <- tb_tumortype$SampleID
      mat_tumortype[is_tcga, 1] <- NA
      colnames(mat_tumortype) <- "CCL"
      ht_ccl <- ComplexHeatmap::Heatmap(mat_tumortype, name = "CCL",
                                        col = col_map_ccl,
                                        na_col = "white",
                                        show_row_names = show_row_names_tumortype,
                                        show_column_names = TRUE,
                                        column_title = NULL,
                                        column_title_rot = 90,
                                        column_title_side = "bottom")

      ht_list <- ht + ht_tcga + ht_ccl
    }
  } else {
    ht_list <- ht
  }

  # category side heatmap
  if (has_category) {
    # category color map
    cat_vector <- df_sample$Category %>%
      unique()
    if (is.null(df_category_color)) {
      if (length(cat_vector) > length(col_vector)) {
        stop('The number of categorys is greater than the number of discrete colors implemented (', length(col_vector), ')')
      }
      col_map_cat <- structure(col_vector[1:length(cat_vector)],
                               names = cat_vector)
    } else {
      if (length(setdiff(c('Category', 'Color'), colnames(df_category_color))) > 0) {
        stop('df_category_color must have columns: Category, Color')
      }

      # convert factors to strings
      if (is.factor(df_category_color$Category)) {
        df_category_color$Category <- as.character(df_category_color$Category)
      }
      if (is.factor(df_category_color$Color)) {
        df_category_color$Color <- as.character(df_category_color$Color)
      }

      if (length(setdiff(cat_vector, df_category_color$Category)) > 0) {
        stop('some Categories are not included in df_category_color')
      }
      col_map_cat <- structure(df_category_color$Color[df_category_color$Category %in% cat_vector],
                               names = df_category_color$Category[df_category_color$Category %in% cat_vector])
    }

    # category side heatmap
    mat_cat <- dplyr::as_tibble(df_sample) %>%
      dplyr::select(Category) %>%
      as.matrix()
    rownames(mat_cat) <- df_sample$SampleID
    cat_name <- ifelse(is.null(category_name),
                       "Category",
                       category_name)
    colnames(mat_cat) <- cat_name
    ht_cat <- ComplexHeatmap::Heatmap(mat_cat, name = cat_name,
                                      col = col_map_cat,
                                      na_col = "white",
                                      show_row_names = show_row_names_cat,
                                      show_column_names = TRUE,
                                      column_title = NULL,
                                      column_title_rot = 90,
                                      column_title_side = "bottom")
    ht_list <- ht_list + ht_cat
  }

  ComplexHeatmap::draw(ht_list, newpage = FALSE)
}
