#' The available methods for statistical analysis.
#' 
#' @return A character vector with the available methods for statistical
#' analysis
stat_analysis_methods <- function() {
  c('limma', 'ttest')
}

#' Differential analysis of protemic data set with different methods.
#' 
#' @param p_df The data to run differential statistical analysis on. 
#'   Of class `proteomics_data`.
#' @param method The method used for differential analysis. See
#'   `stat_analysis_methods()` for available options.
#' @param contrasts The contrasts between condition to test in differential
#'   analysis. Expected to be a data frame with two columns, a and b, where 
#'   each row corresponds to a contrast to test.
#' @param adjust_method The method used for correcting for multiple hypothesis
#'  testing. See `p.adjust.methods` for available options.
#' @return A data frame with calculated log fold changes and adjusted pvalues
#'   for each specified contrast.
#'
#' @export
diff_expr <- function(
    p_df,
    method,
    contrasts,
    adjust_method = 'BH') {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  # extract from input data
  df_wide <- p_df %>%
    tidyr::pivot_wider(names_from = 'label', values_from = 'LFQ', values_fill = NA)
  mx <- df_wide %>%
    dplyr::select(tidyselect::where(is.numeric)) %>%
    as.matrix
  rownames(mx) <- df_wide %>% 
    dplyr::select(!tidyselect::where(is.numeric)) %>% 
    dplyr::pull
  annotation <- attr(p_df, 'annotation')
  # get groups from data and annotation
  cols <- mx %>% colnames
  groups <- annotation$group[order(match(annotation$label, cols))]
  
  if (method == 'limma') {
    # transform contrasts into limma format
    groups <- factor(make.names(groups))
    contrasts <- contrasts %>%
      dplyr::transmute(stringr::str_glue('{make.names(a)} - {make.names(b)}')) %>%
      dplyr::pull
    # run analysis
    out <- .stat_analysis_limma(mx, groups, contrasts)
  } else if (method == 'ttest') {
    out <- .stat_analysis_ttest(mx, group, contrasts)
  }
  
  # adjust pvalues and return output
  out %>%
    dplyr::mutate(pval_adj = stats::p.adjust(pval, adjust_method))
}


#' Differential analysis using LIMMA
#' 
#' @param mx An expression matrix.
#' @param groups A vector indicating which columns belong to which group
#' @param contrasts A vector of contrasts in limma format.
#' @return A data frame with raw pvalues.
.stat_analysis_limma <- function(
    mx,
    groups,
    contrasts) {
  # create design matrix and fix colnames
  design_mx <- stats::model.matrix(~ 0 + groups)
  colnames(design_mx) <- colnames(design_mx) %>%
    stringr::str_replace(deparse(substitute(groups)), '')
  # create contrast matrix
  contrast_mx <- limma::makeContrasts(contrasts = contrasts, levels = design_mx)
  # run limma
  fit <- limma::lmFit(mx, design_mx)
  contr <- limma::contrasts.fit(fit, contrast_mx)
  bayes <- limma::eBayes(contr)
  # parse output
  colnames(contrast_mx) %>%
    purrr::map(function(contr) {
      limma::topTable(bayes, coef=contr, number=Inf) %>%
        tibble::rownames_to_column('id') %>%
        dplyr::mutate(contrast = contr)
    }) %>% 
    dplyr::bind_rows %>%
    dplyr::rename(logfc = logFC,
                  pval = P.Value) %>%
    dplyr::select(contrast, id, logfc, pval)
}


#' Differential analysis using T-tests.
#' 
#' @param mx An expression matrix.
#' @param groups A vector indicating which columns belong to which group
#' @param contrasts The contrasts between condition to test in differential
#'   analysis. Expected to be a data frame with two columns, a and b, where 
#'   each row corresponds to a contrast to test.
#' @return A data frame with raw pvalues.
.stat_analysis_ttest <- function(
    mx, 
    groups, 
    contrasts) {
  # for each contrast
  contrasts %>%
    purrr::pmap(function(a, b) {
      # extract groups from matrix
      mx1 <- mx[,which(groups == a)]
      mx2 <- mx[,which(groups == b)]
      # run tests
      p_values <- 1:nrow(mx1) %>%
        purrr::map(function(n) {
          # check if enough observations
          if (sum(!is.na(mx1[n,])) < 2 & sum(!is.na(mx2[n,])) < 2) {
            return(NA_real_)
          }
          tryCatch({
            stats::t.test(mx1[n,], mx2[n,]) %>%
              broom::glance %>%
              dplyr::pull(p.value)
          }, error = function(e) { return(NA_real_) })
        }) %>%
	    purrr::flatten_dbl
      
      # get logfcs
      logfc <- unname(rowMeans(mx1, na.rm = T) - rowMeans(mx2, na.rm = T))
      # assemble output
      tibble::tibble(contrast = stringr::str_glue('{make.names(a)} - {make.names(b)}'),
                     id = rownames(mx),
                     logfc = logfc,
                     pval = p_values)
    }) %>%
    dplyr::bind_rows
}

#' Construct contrasts against control group.
#'
#' @param p_df Proteomics data.
#' @param cntrl_group character, the control group.
#' @return data frame, the contrasts against the control group in the expected format.
#'
#' @export
construct_contrasts_control <- function(p_df, cntrl_group) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  all_groups <- unique(attr(p_df, 'annotation')$group)
  
  if (!cntrl_group %in% all_groups) {
    stop('Control group not found in protemic data set.',
         call. = FALSE)
  }
  
  tibble::tibble(a = all_groups[!all_groups == cntrl_group],
                 b = cntrl_group)
}

#' Construct all pairwise contrasts.
#'
#' @param p_df Proteomics data.
#' @return data frame, the constrasts in the expected format.
#'
#' @export
construct_contrasts_all <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  all_groups <- unique(attr(p_df, 'annotation')$group)
  
  utils::combn(all_groups, 2) %>% 
    t %>% 
    tibble::as_tibble(.name_repair = 'unique') %>% 
    dplyr::rename(a = '...1', b = '...2')
}

#' Assign comparison type.
#'
#' Construct info whether comparisons are based on value-value, value-imput, 
#' or imput-imput, for each protein and each contrasts given.
#'
#' In order for this comparison to be accurate, the proteomics data from after
#' filtering, but before imputation should be used.
#'
#' @param p_df Proteomics data.
#' @param contrasts data frame, the contrasts for which the type should be determined.
#'   Expected in the same format as for `diff_expr`.
#' @return Data frame with the contrast types.
#'
#' @export
diff_type <- function(p_df, contrasts) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  contrasts %>%
    purrr::pmap(function(a,b) {
      p_df %>%
        filter_data_by_group(c(a, b)) %>%
        dplyr::inner_join(attr(p_df, 'annotation'), by = 'label') %>%
        dplyr::select(id, group) %>%
        dplyr::distinct %>%
        dplyr::group_by(id) %>%
        dplyr::summarise(diff_type = dplyr::case_when(all(c(a, b) %in% group) ~ 'value_value',
                                                      a %in% group ~ 'value_imput',
                                                      b %in% group ~ 'imput_value',
                                                      TRUE ~ NA_character_),
                  contrast = stringr::str_glue('{make.names(a)} - {make.names(b)}'),
                  .groups = 'drop')
      
        
    }) %>%
    dplyr::bind_rows %>%
    tidyr::pivot_wider(names_from = contrast, values_from = diff_type, 
                       values_fill = 'imput_imput') %>%
    tidyr::pivot_longer(-id, names_to = 'contrast', values_to = 'diff_type') %>%
    dplyr::select(contrast, id, diff_type)
}
