#' The available filtering methods.
#'
#' @return A character vector of available filtering methods.
#'
#' @export
filtering_methods <- function() {
  c('none', 'whole_dataset', 'any_group', 'each_group')
}

#' Filtering proteomics data set with different methods
#'
#' TODO WRITE DETAILS
#'
#' @param p_df The data to run filtering on. Of class `proteomics_data`.
#' @param method The method used for filtering. See `filtering_methods()` for
#'   available option.
#' @param threshold The threshold how many samples need to be present in the
#'   subset the different methods look at to retain a data point.
#' @return A filtered object of class `proteomics_data`.
#'
#' @export
filtering <- function(
    p_df,
    method,
    threshold) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  method <- match.arg(method, choices=filtering_methods())
  if (method %in% c('whole_dataset', 'any_group', 'each_group')) {
    if (!all(threshold >=0, threshold <= 1, is.numeric(threshold), length(threshold) == 1)) {
      stop('Threshold must be a numeric value between zero and one.',
           call. = FALSE)
    }
  }
  
  # extract data
  data <- tibble::as_tibble(p_df)
  annotation <- attr(p_df, 'annotation')
  
  # perform computation
  if (method == 'none') {
    data_filtered <- data
  } else if (method == 'whole_dataset') {
    data_filtered <- .filtering_whole_dataset(data, threshold)
  } else if (method == 'any_group') {
    data_filtered <- .filtering_any_group(data, threshold, annotation)
  } else if (method == 'each_group') {
    data_filtered <- .filtering_each_group(data, threshold, annotation)
  }
  
  # check if all groups are still present between data and annotation
  if (!all(annotation$label %in% unique(data_filtered$label))) {
    missing_labels <- annotation$label[!annotation$label %in% unique(data_filtered$label)]
    warning('Filtering removed the following labels completely from the data:',
            stringr::str_glue(' {missing_labels}'))
    annotation <- annotation %>%
      dplyr::filter(!label %in% missing_labels)
  }
  
  # reform to object
  proteomics_data(
    data_filtered, annotation, 
    has_tech_repl = attr(p_df, 'has_tech_repl'),
    is_log2 = TRUE)
}

#' Filtering data set by presence of values across whole data set
#'
#' @param data The data frame to run filtering on.
#' @param threshold The threshold how many samples need to be present in the
#'   subset the different methods look at to retain a data point.
#' @return A filtered data frame.
.filtering_whole_dataset <- function(data, threshold) {
  n_samples <- data %>% dplyr::pull(label) %>% unique() %>% length()
  data %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(present_samples = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(present_samples/n_samples >= threshold)
}

#' Filtering data set by presence of values in at least one group in the data
#' set
#'
#' @param data The data frame to run filtering on.
#' @param threshold The threshold how many samples need to be present in the
#'   subset the different methods look at to retain a data point.
#' @param annotation A data frame that contains information about which samples
#'   in the input data frame belong to which condition.
#' @return A filtered data frame.
.filtering_any_group <- function(data, threshold, annotation) {
  data %>%
    # join with annotation
    dplyr::inner_join(annotation %>%
                        dplyr::group_by(group) %>%
                        dplyr::mutate(n_group = dplyr::n()) %>%
                        dplyr::ungroup(),
                      by='label') %>%
    # calculate ratio of present sample for each peptide for each group
    dplyr::group_by(id, group) %>%
    dplyr::mutate(present_ratio = dplyr::n()/n_group) %>%
    dplyr::ungroup() %>%
    # get best ratio for each peptide
    dplyr::group_by(id) %>%
    dplyr::mutate(best_present_ratio = max(present_ratio)) %>%
    dplyr::ungroup() %>%
    # filter by best ratio
    dplyr::filter(best_present_ratio >= threshold)
}

#' Filtering data set by presence of values for each group in the data
#'
#' @param data The data frame to run filtering on.
#' @param threshold The threshold how many samples need to be present in the
#'   subset the different methods look at to retain a data point.
#' @param annotation A data frame that contains information about which samples
#'   in the input data frame belong to which condition.
#' @return A filtered data frame.
.filtering_each_group <- function(data, threshold, annotation) {
  df <- data %>%
    # join with annotation
    dplyr::inner_join(annotation %>%
					    dplyr::group_by(group) %>%
					    dplyr::mutate(n_group = dplyr::n()) %>%
					    dplyr::ungroup(),
					  by='label') %>%
    # calculate ratio of present sample for each peptide for each group
    dplyr::group_by(id, group) %>%
    dplyr::mutate(present_ratio = dplyr::n()/n_group) %>%
    dplyr::ungroup() %>%
    # filter by ratio
    dplyr::filter(present_ratio >= threshold)
}
