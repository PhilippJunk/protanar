#' The OOP framework for storing data in this analysis pipeline.
#' 
#' Contains constructors for and manipulations of the `proteomics_data` class.

###############################################################################
# class definitions, as suggested as best practice in "Advanced R"

#' Low level constructor for `proteomics_data` class
#' 
#' @param df The data frame containing the proteomics data.
#' @param annotation A valid annotation data frame.
#' @param has_tech_repl A boolean value indicating whether the data has been
#'   log2 transformed.
#' @return A `proteomics_data` object.
new_proteomics_data <- function(df, annotation, has_tech_repl) {
  stopifnot(is.data.frame(df))
  stopifnot(is.data.frame(annotation))
  stopifnot(is.logical(has_tech_repl))
  stopifnot(is.character(df$id), is.character(df$label), 
            is.character(annotation$label), is.character(annotation$group))
  stopifnot(is.numeric(df$LFQ), is.numeric(annotation$biol_repl))
  if (has_tech_repl) {
    stopifnot(is.numeric(annotation$tech_repl))
  }

  proteomics_data <- structure(
    df,
    annotation = annotation,
    class = c('proteomics_data', 'tbl_df', 'tbl', 'data.frame'),
    has_tech_repl = has_tech_repl)

  return(proteomics_data)
}

#' Validator for `proteomics_data` class
#'
#' @param proteomics_data A `proteomics_data` object.
#' @return A validated `proteomics_data` object.
validate_proteomics_data <- function(proteomics_data) {
  df <- tibble::as_tibble(proteomics_data)
  annotation <- attr(proteomics_data, 'annotation')
  has_tech_repl <- attr(proteomics_data, 'has_tech_repl')

  if (!inherits(proteomics_data, 'proteomics_data')) {
    stop('Please provide an object of class `proteomics_data` as input. ', 
         'For more information, see `?proteomics_data`.',
         call. = FALSE)
  }
  
  if (!is.data.frame(annotation)) {
    stop(
      '`annotation` must be a `data.frame`.',
      call. = FALSE
    )
  }

  if (!is.logical(has_tech_repl)) {
    stop(
      '`has_tech_repl` must be a logical vector.',
      call. = FALSE
    )
  }
  
  if (!all(colnames(df) == c('id', 'label', 'LFQ'))) {
    stop(
      'Data must consists of the columns `id`, `label` and `LFQ`.',
      call. = FALSE
    )
  }
  
  if (has_tech_repl) {
    if (!all(colnames(annotation) == c('label', 'group', 'biol_repl', 'tech_repl'))) {
      stop('Annotation is expected to have the columns `label`, `group`, `biol_repl` and `tech_repl`.', 
           call. = FALSE)
    }
  } else {
    if (!all(colnames(annotation) == c('label', 'group', 'biol_repl'))) {
      stop('Annotation is expected to have the columns `label`, `group` and `biol_repl`.', 
           call. = FALSE)
    }
  }
  
  if (any(is.na(df$LFQ))) {
    stop('LFQ intensities must be non-missing.',
         call. = FALSE)
    
  }
      
  if (any(df$LFQ < 0)) {
    stop('LFQ intensities must be non-missing and greater than zero.',
         call. = FALSE)
  }
  
  if (!all(unique(df$label) %in% annotation$label)) {
    stop('Not all `label`s in data are present in `annotation`.', 
         call. = FALSE)
  }
  
  if (!all(annotation$label %in% unique(df$label))) {
    stop('Not all `label`s in `annotation` are present in data', 
         call. = FALSE)
  }
  
  if (any(duplicated(annotation$label))) {
    stop('`label`s in `annotation` must be unique.', 
         call. = FALSE)
  }
  
  return(proteomics_data)
}

#' Construct a proteomics data object.
#' 
#' @param df The data frame with proteomics data.
#' @param annotation The data frame with annotation for the proteomics data.
#' @param has_tech_repl logical, whether the data and annotation has technical replicates.
#' @param is_log2 logical, whether the proteomics data is already log2 transformed.
#' @param df_id character, the column name in `df` that uniquely identifies different proteins (default is id).
#' @param df_label character, the column name in `df` that assignes a unique label to each sample (default is label).
#' @param df_LFQ character, the column name in `df` for quantitative values (defaults is LFQ).
#' @param annotation_label character, the column name in `annotation` that contains the sample labels used in `df$label` (default is label).
#' @param annotation_group character, the column name in `annotation` that contains the groups the each sample belongs to (default is group). Multiple samples can belong to the same group.
#' @param annotation_biol_repl character, the column name in `annotation` that contains the identifier of the biological replicate of the sample (default is biol_repl).
#' @param annotation_tech_repl character, the column name in `annotation` that contains the identifier of the technical replicate of the sample (default is tech_repl).
#' @return A `proteomics_data` object.
#'
#' @export
proteomics_data <- function(
    df,
	annotation,
	has_tech_repl,
	is_log2,
    df_id = 'id',
	df_label = 'label',
	df_LFQ = 'LFQ', 
    annotation_label = 'label',
	annotation_group = 'group', 
    annotation_biol_repl = 'biol_repl',
	annotation_tech_repl ='tech_repl'
) {
  
  df <- df %>% 
    dplyr::select({df_id}, {df_label}, {df_LFQ}) %>%
    dplyr::rename(id = {df_id},
				  label = {df_label},
				  LFQ = {df_LFQ}) %>%
    dplyr::select(id, label, LFQ) %>%
    dplyr::mutate(id = make.names(id),
				  label = make.names(label))
  if (!is_log2) {
    stopifnot(all(df$LFQ > 2))
    df <- df %>%
      dplyr::mutate(LFQ = log2(LFQ))
  }
  if (has_tech_repl) {
    annotation <- annotation %>%
      dplyr::select({annotation_label}, {annotation_group}, 
                    {annotation_biol_repl}, {annotation_tech_repl}) %>%
	  dplyr::rename(label = {annotation_label},
                    group = {annotation_group},
                    biol_repl = {annotation_biol_repl},
                    tech_repl = {annotation_tech_repl}) %>%
      dplyr::select(label, group, biol_repl, tech_repl)
  } else {
    annotation <- annotation %>%
      dplyr::select({annotation_label}, {annotation_group}, 
					{annotation_biol_repl}) %>%
      dplyr::rename(label = {annotation_label},
					group = {annotation_group},
					biol_repl = {annotation_biol_repl}) %>%
      dplyr::select(label, group, biol_repl)
  }
  annotation %>% 
    dplyr::mutate(group = make.names(group),
				  label = make.names(label))
  
  validate_proteomics_data(new_proteomics_data(
    df, annotation,has_tech_repl))
}

###############################################################################
# helper functions

#' Join proteomics data with associated annotation.
#'
#' @param p_df Proteomics data object.
#' @return Data frame of proteomics data joined with annotation.
#'
#' @export
join_annotation <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  
  p_df %>%
    dplyr::inner_join(attr(p_df, 'annotation'), by = 'label')
}

#' Collapse technical replicates in proteomics data.
#'
#' @param p_df Proteomics data with technical replicates.
#' @return Proteomics data.
#'
#' @export
collapse_tech_repl <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  if (!attr(p_df, 'has_tech_repl')) {
    stop('The `proteomics_data` object must have technical replicates in order to collapse them.',
         call. = FALSE)
  }
  
  annotation <- attr(p_df, 'annotation') %>%
    dplyr::mutate(new_label = stringr::str_c(group, biol_repl, sep='_'))
  
  p_df <- p_df %>%
    # transform back from log2 for collapsing
    dplyr::mutate(LFQ = 2^LFQ) %>%
    dplyr::inner_join(annotation, by = 'label') %>%
    dplyr::group_by(id, new_label) %>%
    dplyr::summarise(LFQ = stats::median(LFQ), .groups='drop')
  
  annotation <- annotation %>%
    dplyr::select(new_label, group, biol_repl) %>%
    dplyr::distinct
    
  proteomics_data(p_df, annotation, has_tech_repl = FALSE, is_log2 = FALSE,
                  df_label = 'new_label', annotation_label = 'new_label')
} 

# filter data by group and id
# df_filter expected to have group and id column!

#' Filer proteomics data by groups and ids.
#'
#' This filtering can be useful to filter a data set against the output of a
#' differential analysis, for example when removing proteins from a data set that
#' have also been identified in a negative control.
#' 
#' @param p_df Proteomics data.
#' @param df_filter Dataframe, expected to have two columns called group and id.
#' @return Filtered proteomics data.
#'
#' @export
filter_data_by_group_and_id <- function(p_df, df_filter) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  df_filter <- df_filter %>%
    dplyr::select(id, group)

  # filter data
  p_df <- p_df %>%
    dplyr::inner_join(attr(p_df, 'annotation'), by='label') %>%
    dplyr::inner_join(df_filter, by = c('id', 'group'))
  # filter annotation
  annotation <- attr(p_df, 'annotation') %>%
    dplyr::inner_join(df_filter %>% dplyr::select(group) %>% dplyr::distinct,
					  by = 'group') 
  
  proteomics_data(
    p_df, annotation, 
    has_tech_repl = attr(p_df, 'has_tech_repl'), is_log2 = TRUE)
}

#' Filter proteomics data by groups.
#'
#' @param p_df Proteomics data.
#' @param groups Vector of groups that the data set will be filtered by.
#' @param keep Logical, whether members of groups should be kept in the
#'   data set or discarded (default is TRUE, so specified groups will be kept).
#' @return Filtered proteomics data.
#'
#' @export
filter_data_by_group <- function(p_df, groups, keep=TRUE) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  annotation <- attr(p_df, 'annotation')

  # reverse group selection of keep=FALSE
  if (!keep) {
	all_groups <- annotation$group %>% unique
	groups <- all_groups[!all_groups %in% groups]
  }

  # filter data
  p_df <- p_df %>%
    dplyr::inner_join(attr(p_df, 'annotation'), by='label') %>%
    dplyr::filter(group %in% groups)
  # filter annotation
  annotation <- annotation %>%
    dplyr::filter(group %in% groups)
  
  proteomics_data(
    p_df, annotation, 
    has_tech_repl = attr(p_df, 'has_tech_repl'), is_log2 = TRUE)
}

#' Filter proteomics data by samples.
#'
#' @param p_df Proteomics data.
#' @param samples Vector of samples that the data set will be filered by.
#' @param keep logical, whether members of samples should be kept in the
#'   data set or discarded (default is TRUE, so specified samples will be kept).
#' @return Filtered proteomics data.
#'
#' @export
filter_data_by_samples <- function(p_df, samples, keep=TRUE) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  annotation <- attr(p_df, 'annotation')

  # reverse sample selection if keep=FALSE
  if (!keep) {
	all_samples <- annotation$labels %>% unique
	samples <- all_samples[!all_samples %in% samples]
  }
  
  # filter data
  p_df <- p_df %>%
    dplyr::filter(!label %in% samples)
  # filter annotation
  annotation <- annotation %>%
    dplyr::filter(!label %in% samples)
  
  proteomics_data(
    p_df, annotation,
    has_tech_repl = attr(p_df, 'has_tech_repl'), is_log2 = TRUE)
}

#' Deconstruct groups in proteomics data.
#'
#' TODO WHAT DOES THIS DO AGAIN?
#'
#' @param p_df Proteomics data.
#' @return TODO
#'
#' @export
deconstruct_groups <- function(p_df) {
  # check inputs
  p_df <- validate_proteomics_data(p_df)
  annotation <- attr(p_df, 'annotation')

  sets <- p_df %>%
    dplyr::inner_join(annotation, by='label') %>%
    dplyr::mutate(group = factor(group) %>% forcats::fct_infreq %>% forcats::fct_rev) %>%
    dplyr::select(id, group) %>%
    dplyr::distinct %>%
    dplyr::arrange(group) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(set = stringr::str_c(group, collapse = '__'),
					 .groups = 'drop')
  
  sets %>%
    dplyr::count(set) %>% 
    dplyr::arrange(-n) %>%
    dplyr::pull(set) %>%
    rlang::set_names %>%
    purrr::map(function(s) {
      sets %>%
        dplyr::filter(set == s) %>%
        dplyr::pull(id)
    })
}
