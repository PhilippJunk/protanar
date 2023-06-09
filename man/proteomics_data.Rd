% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_obj.R
\name{proteomics_data}
\alias{proteomics_data}
\title{Construct a proteomics data object.}
\usage{
proteomics_data(
  df,
  annotation,
  has_tech_repl,
  is_log2,
  df_id = "id",
  df_label = "label",
  df_LFQ = "LFQ",
  annotation_label = "label",
  annotation_group = "group",
  annotation_biol_repl = "biol_repl",
  annotation_tech_repl = "tech_repl"
)
}
\arguments{
\item{df}{The data frame with proteomics data.}

\item{annotation}{The data frame with annotation for the proteomics data.}

\item{has_tech_repl}{logical, whether the data and annotation has technical replicates.}

\item{is_log2}{logical, whether the proteomics data is already log2 transformed.}

\item{df_id}{character, the column name in \code{df} that uniquely identifies different proteins (default is id).}

\item{df_label}{character, the column name in \code{df} that assignes a unique label to each sample (default is label).}

\item{df_LFQ}{character, the column name in \code{df} for quantitative values (defaults is LFQ).}

\item{annotation_label}{character, the column name in \code{annotation} that contains the sample labels used in \code{df$label} (default is label).}

\item{annotation_group}{character, the column name in \code{annotation} that contains the groups the each sample belongs to (default is group). Multiple samples can belong to the same group.}

\item{annotation_biol_repl}{character, the column name in \code{annotation} that contains the identifier of the biological replicate of the sample (default is biol_repl).}

\item{annotation_tech_repl}{character, the column name in \code{annotation} that contains the identifier of the technical replicate of the sample (default is tech_repl).}
}
\value{
A \code{proteomics_data} object.
}
\description{
Construct a proteomics data object.
}
