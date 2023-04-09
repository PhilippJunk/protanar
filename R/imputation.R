#' Imputation methods
#'
#' @description
#' `imputation_methods()` returns the names of all available imputation methods.
#'
#' `imputation_methods_mixed` returns the names of all available imputation
#' methods that use mixed imputation. 
#'
#' When a mixed imputation strategy is used, the methods returned by
#' `imputation_methods_mar` and `imputation_methods_mnar` are available for
#' "missing at random" (MAR) imputation and "missing not at random" (MNAR)
#' imputation, respectively.
#'
#' @export
imputation_methods <- function() {
  c(imputation_methods_mar(), imputation_methods_mnar(),
    imputation_methods_mixed())
}

#' @rdname imputation_methods
#' @export
imputation_methods_mixed <- function() {
  c('mixed_row', 'mixed_sample')
}

#' @rdname imputation_methods
#' @export
imputation_methods_mar <- function() {
  c('MLE', 'bpca', 'knn')
}

#' @rdname imputation_methods
#' @export
imputation_methods_mnar <- function() {
  c('QRILC', 'MinDet', 'MinProb', 'zero', 'min')
}

#' Impute expression data set with different methods.
#'
#' @param p_df The data to run filtering on. Of class `proteomics_data`.
#' @param method The method used for imputation. See `imputation_methods()` for
#'   available options.
#' @param mar_method The method used for the imputation of missing-at-random
#'   (MAR) missing value if a mixed imputation approach has been chosen. See
#'   `imputation_methods_mar()` for available options.
#' @param mnar_method The method used for the imputation of missing-not-at-
#'   random (MNAR) missing values if a mixed imputation approach has been
#'   chosen. See `imputation_methods_mnar()` for available options.
#' @param ... Futher arguments passed down to helper functions.
#' @return An object of class `proteomics_data` with imputed values.
#'
#' @export
imputation <- function(
    p_df,
    method,
    mar_method,
    mnar_method,
    ...
) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  method = match.arg(method, choices = imputation_methods())
  if (method %in% imputation_methods_mixed()) {
    mar_method <- match.arg(mar_method, choices = imputation_methods_mar())
    mnar_method <- match.arg(mnar_method, choices = imputation_methods_mnar())
  }

  # extract from input data
  df_wide <- p_df %>%
    tidyr::pivot_wider(names_from = 'label', values_from = 'LFQ', values_fill = NA)
  mx <- df_wide %>%
    dplyr::select(tidyselect::where(is.numeric)) %>%
    as.matrix
  cols <- mx %>% colnames
  annotation <- attr(p_df, 'annotation')

  # perform imputation
  if (! method %in% c('mixed_row', 'mixed_sample')) {
    mx_imput <- MsCoreUtils::impute_matrix(mx, method, ...)
  } else {
    # obtain groups from data
    groups <- annotation$group[order(match(annotation$label, cols))]
    if (method == 'mixed_row') {
      # calculate MAR vector
      mar_vct <- .get_mar_vector(mx, groups)
      # perform imputation
      mx_imput <- MsCoreUtils::impute_mixed(mx, randna = mar_vct, mar = mar_method,
                                            mnar = mnar_method, ...)
    } else if (method == 'mixed_sample') {
      # calculate MAR matrix
      mar_mx <- .get_mar_matrix(mx, groups)
      # perform imputation
      mx_imput <- .imputation_mixed_sample(mx, mar_mx, mar_method, mnar_method, ...)
    }
  }
  if (any(is.na(mx_imput))) {
    stop('Imputation unsuccessful. There are still missing values in data frame. ',
         'Please try a different method, depending on the data, some methods are not applicable.',
         call. = FALSE)
  }

  # reassemble data frame
  data_imputed <- dplyr::bind_cols(
    df_wide %>% dplyr::select(!tidyselect::where(is.numeric)),
    tibble::as_tibble(mx_imput)
  ) %>% 
    tidyr::pivot_longer(cols = tidyselect::all_of(cols),
						names_to = 'label', values_to = 'LFQ')
  
  # reform to object
  proteomics_data(
    data_imputed, annotation, 
    has_tech_repl = attr(p_df, 'has_tech_repl'),
    is_log2 = TRUE)
}


#' Sample-wise mixed imputation
#'
#' Impute expression matrix with a mixed imputation using a sample-wise
#' separation in MAR and MNAR missing values
#'
#' @param mx A expression matrix with missing values.
#' @param mar_mx A logical matrix of the same dimensions as mx. TRUE values
#'   are considered MAR.
#' @param mar_method The method used for imputing MAR missing values.
#' @param mnar_method The method used for imputing MNAR missing values.
#' @param ... Other parameters passed on to mar_method and mnar_method.
#' @return A matrix of the same dimensions as the input with imputed values.
.imputation_mixed_sample <- function(
    mx, mar_mx, mar_method, mnar_method, ...) {
  mx <- MsCoreUtils::impute_matrix(mx, mnar_method, ...)
  mx[mar_mx] <- NA
  mx <- MsCoreUtils::impute_matrix(mx, mar_method, ...)
  mx
}


#' Obtain matrix of MAR values.
#'
#' Missing at random values are defined as observations where not all values
#' for a specific group/condition are missing, but only some.
#'
#' @param mx A expression matrix with missing values.
#' @param groups A vector indicating which columns belong to the same group.
#' @return A logical matrix of the same dimension as mx. TRUE values are
#'   are considered MAR.
.get_mar_matrix <- function(mx, groups) {
  out <- matrix(F, nrow=nrow(mx), ncol=ncol(mx))

  # iterate groups and rows
  for (gr in unique(groups)) {
    g <- groups == gr
    for (r in 1:nrow(mx)) {
      if (any(is.na(mx[r,g])) & !all(is.na(mx[r,g]))) {
        out[r,g] <- is.na(mx[r,g])
      }
    }
  }

  return(out)
}


#' Obtain vector of MAR rows.
#'
#' Missing at random values are defined as observations where not all values
#' for a specific group/condition are missing, but only some.
#'
#' @param mx A expression matrix with missing values.
#' @param groups A vector indicating which columns belong to the same group.
#' @return A logical of length nrow(mx). TRUE values are considered MAR.
.get_mar_vector <- function(mx, groups) {
  out <- !vector(mode = 'logical', length = nrow(mx))
  # iterate rows
  for (r in 1:nrow(mx)) {
    # check if any group is MNAR
    any_mnar <- unique(groups) %>%
      purrr::map(function(x) {
        g <- groups == x
        all(is.na(mx[r,g]))
      }) %>%
	  purrr::flatten_lgl %>%
      any
    if (any_mnar) {
      out[r] <- F
    }
  }
  return(out)
}
