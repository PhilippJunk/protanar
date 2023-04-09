# ideas:
# Volcano plots
# heatmap?
# UMAP?

#' PCA plot for proteomics data.
#'
#' @param p_df Proteomics data.
#'
#' @export
vis_pca <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  df_wide <- p_df %>%
    tidyr::pivot_wider(names_from = id, values_from = LFQ, values_fill = 0)
  
  pca <- df_wide %>%
    dplyr::select(tidyselect::where(is.numeric)) %>%
    dplyr::select(tidyselect::where(~ sd(.x) > 0)) %>%
    stats::prcomp(scale=T, center=T)
  
  pc1_varexpl <- round(summary(pca)$importance[2,1] * 100, 2)
  pc2_varexpl <- round(summary(pca)$importance[2,2] * 100, 2)
  
  pca %>%
    broom::augment(df_wide) %>%
    dplyr::inner_join(attr(p_df, 'annotation'), by= 'label') %>%
    ggplot2::ggplot(ggplot2::aes(x=.fittedPC1, y=.fittedPC2, color=group)) +
      ggplot2::geom_point() +
      ggplot2::labs(x = stringr::str_c('PC1 (', pc1_varexpl, '%)'),
           y  = stringr::str_c('PC2 (', pc2_varexpl, '%)')) +
      ggplot2::theme_bw() +
      NULL
}


#' Scatter plot for proteomics data quality control.
#'
#' Plots scatter plots of all pairwise comparisons between two samples in the data set.
#'
#' @param p_df Proteomics data.
#'
#' @export
vis_qc_scatter <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  labels <- attr(p_df, 'annotation') %>%
    dplyr::pull(label) %>%
    unique %>%
    sort
  min_LFQ <- min(p_df$LFQ)
  max_LFQ <- max(p_df$LFQ)
  
  tidyr::expand_grid(l1 = labels, l2 = labels) %>%
    purrr::pmap(function(l1, l2) {
      df1 <- p_df %>% 
        dplyr::filter(label == l1) %>%
        dplyr::rename(LFQ_1 = LFQ)
      df2 <- p_df %>%
        dplyr::filter(label == l2) %>%
        dplyr::rename(LFQ_2 = LFQ)
      df <- dplyr::inner_join(df1, df2, by='id') 
      corr_LFQ <- stats::cor(df$LFQ_1, df$LFQ_2) %>%
        round(2)
      df %>%
        ggplot2::ggplot(ggplot2::aes(x = LFQ_2, y = LFQ_1)) +
          ggplot2::geom_point() +
          ggplot2::scale_x_continuous(limits = c(min_LFQ, max_LFQ)) +
          ggplot2::scale_y_continuous(limits = c(min_LFQ, max_LFQ)) +
          ggplot2::geom_abline(slope=1, intercept=0) +
          ggplot2::labs(title = stringr::str_glue('cor: {corr_LFQ}'), x = l2, y = l1) +
		  NULL
    }) %>%
	patchwork::wrap_plots(ncol = length(labels), nrow = length(labels), byrow=T)
}

#' Histograms for quality control of proteomics data.
#'
#' Plot the distributions of the data for each sample as a histogram.
#' 
#' @param p_df Proteomics data.
#'
#' @export
vis_qc_histo <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  ggplot2::ggplot(p_df, ggplot2::aes(x=LFQ, fill=label)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(label ~ .) +
    ggplot2::theme(legend.position = 'none') +
	NULL
}

#' Plot counts for quality control of proteomics data.
#'
#' Plot the counts of proteins for each sample as a bar chart.
#' 
#' @param p_df Proteomics data.
#'
#' @export
vis_qc_count <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  p_df %>%
    dplyr::count(label, name='count') %>% 
    ggplot2::ggplot(ggplot2::aes(y = label, x = count)) +
      ggplot2::geom_bar(stat='identity', color='black') +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
	  NULL
}

#' Plot upset plot of proteomics data.
#'
#' Visualize which proteins can be found in which groups using an upset plot.
#' 
#' @param p_df Proteomics data.
#'
#' @export
vis_upset <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  annotation <- attr(p_df, 'annotation')
  
  nsets <- annotation$group %>% unique %>% length
  nintersects = 2^nsets
  
  p_df %>%
    dplyr::inner_join(annotation, 'label') %>%
    dplyr::select(group, id) %>%
    dplyr::distinct %>%
    dplyr::mutate(present = 1) %>%
    tidyr::pivot_wider(names_from = group, values_from = present, values_fill = 0) %>% 
    as.data.frame %>% 
    UpSetR::upset(sets = annotation %>% dplyr::pull(group) %>% unique, 
          order.by = "freq",
          nsets = nsets,
          nintersects = nintersects)
}
