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
