#' Artificial example data sets
#' 
#' @description
#' `example_protanar_1` is an artifical example data set that contains
#' data for 3 proteins in two different groups, with each 6 replicates.
#' For each group, there are 3 biological replicates at 2 technical replicates each.
#' There are no missing values in this data set, so it can be used as an example for
#' differential analysis directly.
#'
#' `example_protanar_2` is an artifical example data set that contains
#' data for 3 proteins in two different groups, with each 6 replicates.
#' For each group, there are 3 biological replicates at 2 technical replicates each.
#' There are some missing values in this data set, so it can be used as an example for
#' quality control, data filtering and data imputation.
#'
#' `example_protanar_anno` is the annotation data set that accompanies
#' `example_protanar_1` and `example_protanar_2`.
#'
#' @format ## `example_protanar_1`
#' A data frame with 36 rows and 3 columns:
#' \describe{
#'   \item{id}{Unique ID for each protein.}
#'   \item{label}{Unique label for each sample.}
#'   \item{LFQ}{LFQ value (log2 transformed).}
#' }
"example_protanar_1"

#' @rdname example_protanar_1
#' @format ## `example_protanar_2`
#' A data frame with 30 rows and 3 columns:
#' \describe{
#'   \item{id}{Unique ID for each protein.}
#'   \item{label}{Unique label for each sample.}
#'   \item{LFQ}{LFQ value (log2 transformed).}
#' }
"example_protanar_2"

#' @rdname example_protanar_1
#' @format ## `example_protanar_anno`
#' A data frame with 12 rows and 4 columns:
#' \describe{
#'   \item{label}{Unique label for each sample.}
#'   \item{group}{Unique label for each group.}
#'   \item{biol_repl}{Unique number for each biological replicate within a group.}
#'   \item{tech_repl}{Unique number for each technical replicate within a biological replicate for a specific group.}
#' }
"example_protanar_anno"
