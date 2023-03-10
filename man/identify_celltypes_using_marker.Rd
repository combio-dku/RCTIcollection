% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/celltype_identifiers_collection.R
\name{identify_celltypes_using_marker}
\alias{identify_celltypes_using_marker}
\title{identify_celltypes_using_marker}
\usage{
identify_celltypes_using_marker(
  df_data,
  marker,
  method = NULL,
  N_pca = 15,
  Clustering_resolution = 1,
  N_features = 2000
)
}
\arguments{
\item{df_data}{input cell-by-gene matrix (data frame), column name must be in Hugo symbol}

\item{marker}{it can be either a string or a data frame (marker-matrix). For 'HiCAT', it must be a string representing the full path to the HiCAT marker file. For others, it must be a data frame obtained by running 'get_mkr_mat_from_hicat_mkr_file'.}

\item{method}{method to use. It can be one of 'HiCAT', 'SCINA', 'scSorter', 'scType', 'scCatch', or 'Garnett'.}
}
\value{
a vector of predicted celltypes. Each element corresponds to a row of df_data. For HiCAT, it is a data frame with 3 columns, each corresponds to major-type, minor-type and subset.
}
\description{
This function runs one of the marker-based celltype identifiers, specified by method.
}
