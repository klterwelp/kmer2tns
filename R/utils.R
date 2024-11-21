#' Log Messages with Timestamps
#'
#' This function logs messages with a timestamp.
#'
#' @param ... Anything to be concatenated into the log message
#'
#' @return None
#' @export
log_message <- function(...) {
  msg <- paste(...)
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
}

#' Validate kmer_count_table format
#'
#' @param kmer_counts data.table with two columns: "kmer" and "count".
#'
#' @return TRUE if kmer_count is in the correct format. FALSE if not.
#' @export
#'
#' @examples
#' validate_kmer_count_table(ecoliCounts)
validate_kmer_count_table <- function(kmer_counts) {
  # Check if kmer_counts is a data.table
  if (!data.table::is.data.table(kmer_counts)) {
    stop("kmer_counts must be a data.table.")
  }

  # Check if it has exactly two columns named "kmer" and "count"
  if (!all(c("kmer", "count", "positions") %in% colnames(kmer_counts))) {
    stop("kmer_counts must have three columns: 'kmer', 'count', and 'positions'.")
  }

  # Check if "kmer" values are all characters
  if (!is.character(kmer_counts$kmer)) {
    stop("The 'kmer' column must contain character values.")
  }

  # Check if "count" values are all numeric
  if (!is.numeric(kmer_counts$count)) {
    stop("The 'count' column must contain numeric values.")
  }

  # Check if "positions" values are a list
  if (!is.list(kmer_counts$positions)) {
    stop("The 'count' column must contain numeric values.")
  }
  # If all checks pass
  return(TRUE)
}
