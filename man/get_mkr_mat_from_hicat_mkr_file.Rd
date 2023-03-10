% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/celltype_identifiers_collection.R
\name{get_mkr_mat_from_hicat_mkr_file}
\alias{get_mkr_mat_from_hicat_mkr_file}
\title{Get marker-matrix from HiCAT marker file}
\usage{
get_mkr_mat_from_hicat_mkr_file(mkr_file, taxo_level = "major", tissues = NULL)
}
\arguments{
\item{mkr_file}{full path to the HiCAT marker file.}

\item{taxo_level}{taxonomy level. It can be one of 'major', 'minor' or 'subset'.}

\item{tissues}{a vector of tissues (groups of celltypes) to be exported. Use tissue names in the HiCAT marker file.}
}
\value{
marker-matrix in data frame
}
\description{
This function uses HiCAT marker file to generate marker-matrix to be used for
celltype identifiers other than HiCAT. Marker-matrix is a data frame with its
rows being celltypes and columns being markers. Each entry is 0 or 1 indicating
the inclusion of the marker in a celltype.
}
