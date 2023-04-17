#' Collapse UNIPROT IDs
#'
#' This function takes a vector of UNIPROT IDs as produced by Maxquant
#' (in the column "Majority Protein IDs"), collapses isoforms. If after that,
#' multiple IDs are found, there are different modes to select only one.
#'
#' @param uniprot_ids A vector of UNIPROT IDs
#' @param multiple_mode Determines how to deal with entries that can not be
#'   broken down into a single entry. Possible values are "none", "manual",
#'   and "first" (default is "first").
#' @param manual_file In multiple_mode "manual", the file used for manual
#'   replacements. Required data frame needs to be saved in in a tab separated
#'   file with the columns "input" and "output".
#' @param separator_in The separator between multiple UNIPROT IDs and isoforms
#'   in an input (default ";").
#' @param separator_out In multiple_mode "none", the separator between
#'   multiple UNIPROT IDs when writing output (default ";").
#' @return Vector of collapes UNIPROT IDs
#'
#' @export
collapse_uniprot_ids <- function(uniprot_ids,
                                 multiple_mode = "first",
                                 manual_file,
                                 separator_in = ";",
                                 separator_out = ";") {
  # check inputs
  if (!multiple_mode %in% c("none", "manual", "first")) {
    stop("Given mode not valid: ", multiple_mode)
  }

  if (multiple_mode == "manual") {
    # import data set
    manual <- utils::read.table(manual_file, header = TRUE, sep = "\t")
  }

  # only iterate over unique elements
  uniprot_ids_unique <- unique(uniprot_ids)
  out_unique <- vector("character", length(uniprot_ids_unique))

  for (i in seq_along(uniprot_ids_unique)) {
    id <- uniprot_ids_unique[i]
    out_id <- c()
    # split by input separator
    for (split_id in unlist(stringr::str_split(id, separator_in))) {
      # if starts with CON__, skip
      if (stringr::str_detect(split_id, "CON__")) {
        next()
      }
      # if starts with REV__, skip
      else if (stringr::str_detect(split_id, "REV__")) {
        next()
      }
      # if dash found, extract everything in front of it
      else if (stringr::str_detect(split_id, "-\\d+")) {
        out_id <- c(out_id, stringr::str_extract(split_id, "^[^-]+(?=-\\d+)"))
      }
      # if none of those, extract everything
      else {
        out_id <- c(out_id, split_id)
      }
    }
    # take unique ids
    out_id <- unique(out_id)

    # if more than one UNIPROT id was found, proceed
    # depending on mode
    if (multiple_mode == "none") {
      out_id <- stringr::str_c(out_id, collapse = separator_out)
    } else if (multiple_mode == "first") {
      out_id <- out_id[1]
    } else if (multiple_mode == "manual") {
      out_id_concat <- stringr::str_c(out_id, collapse = separator_out)
      if (length(out_id) == 1) {
        # pass
        invisible()
      } else if (out_id_concat %in% manual$input) {
        out_id <- manual %>%
          dplyr::filter(input == out_id_concat) %>%
          dplyr::pull(output)
      } else {
        warning(
          out_id, " not found in manual replacement. ",
          "Selecting first element."
        )
        out_id <- out_id[1]
      }
    }
    out_unique[i] <- out_id
  }

  # apply unique mapping to input
  dplyr::left_join(data.frame(input = uniprot_ids),
    data.frame(
      input = uniprot_ids_unique,
      output = out_unique
    ),
    by = "input"
  ) %>%
    dplyr::pull(output)
}


#' Map UNIPROT IDs to HGNC gene names
#'
#' Given a vector of UNIPROT IDs , extract corresponding HGNC gene
#' names where applicable.
#'
#' @param uniprot_ids A vector of UNIPROT IDs.
#' @param hgnc_db A data.frame of the HGNC data base. For this function,
#'   expected to contain the columns "hgnc_id", "symbol", and "uniprot_ids".
#' @param multiple_mode Determines how cases are treated where one UNIPROT ID
#'   corresponds to multiple HGNC gene names. Possible values are "first" and
#'   "manual" and "none" (default "first").
#' @param manual_file In multiple_mode "manual", the file used for manual
#'   replacements. Required data frame needs to be saved in in a tab separated
#'   file with the columns "input" and "output".
#' @param separator_in The separator used to split multiple UNIPROT ID entries
#'   per input. Default is ";".
#' @param separator_out The separator used for separating multiple HGNC entries
#'   associated with a single UNIPROT ID, if multiple_mode is "none". The
#'   default is ";".
#' @return Vector of HGNC gene names
#'
#' @export
map_uniprot_hgnc <- function(uniprot_ids,
                             hgnc_db,
                             multiple_mode = "first",
                             manual_file,
                             separator_in = ";",
                             separator_out = ";") {
  # check inputs
  if (!multiple_mode %in% c("none", "manual", "first")) {
    stop("Given mode not valid: ", multiple_mode)
  }

  if (multiple_mode == "manual") {
    # import data set
    manual <- utils::read.table(manual_file, header = TRUE, sep = "\t")
  }

  stopifnot(is.data.frame(hgnc_db))
  hgnc_db <- hgnc_db %>%
    dplyr::select(hgnc_id, symbol, uniprot_ids) %>%
    tidyr::separate_rows(uniprot_ids) %>%
    dplyr::filter(!is.na(uniprot_ids))

  # For element in uniprot_id, find assignment in HGNC
  # There can be multiple HGNC IDs associated with one UNIPROT ID
  # There can also be multiple UNIPROT IDs associated with one HGNC ID

  # make input unique
  uniprot_ids_unique <- unique(uniprot_ids)
  out_unique <- vector("character", length(uniprot_ids_unique))

  # perform mapping
  for (i in seq_along(uniprot_ids_unique)) {
    id <- uniprot_ids_unique[i]
    hits <- c()
    # split by semi-colon
    for (split_id in unlist(stringr::str_split(id, separator_in))) {
      hits <- c(hits, hgnc_db %>%
        dplyr::filter(uniprot_ids == split_id) %>%
        dplyr::pull(symbol))
    }
    # take unique ids
    hits <- unique(hits)

    # if more than one UNIPROT id was found, proceed
    # depending on mode
    if (multiple_mode == "none") {
      hits <- stringr::str_c(hits, collapse = separator_out)
    } else if (multiple_mode == "first") {
      hits <- hits[1]
    } else if (multiple_mode == "manual") {
      hits_concat <- stringr::str_c(hits, collapse = separator_out)
      if (length(hits) == 1) {
        # pass
        invisible()
      } else if (hits_concat %in% manual$input) {
        hits <- manual %>%
          dplyr::filter(input == hits_concat) %>%
          dplyr::pull(output)
      } else {
        warning(
          out_id, " not found in manual replacement. ",
          "Selecting first element."
        )
        hits <- out_id[1]
      }
    }
    out_unique[i] <- hits
  }

  # apply unique mapping to input
  dplyr::left_join(data.frame(input = uniprot_ids),
    data.frame(
      input = uniprot_ids_unique,
      output = out_unique
    ),
    by = "input"
  ) %>%
    dplyr::pull(output)
}
