#' PhiX k-mer counts data
#'
#' Counts of k-mers appearing more than once in the PhiX genome
#'
#' @format
#' A data frame with 401 rows and 2 columns:
#' \describe{
#'   \item{kmer}{8-mer sequence}
#'   \item{count}{number of times seen in the sequence}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/nuccore/9626372>
"PhiXcounts"

#' Ecoli k-mer counts data
#'
#' Counts of k-mers appearing more than once in the Ecoli genome
#'
#' @format
#' A data frame with 999708 rows and 2 columns:
#' \describe{
#'   \item{kmer}{12-mer sequence}
#'   \item{count}{number of times seen in the sequence}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/nuccore/AP022815.1>
"ecoliCounts"
