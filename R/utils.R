#' Log Messages with Timestamps
#'
#' This function logs messages with a timestamp.
#'
#' @param msg A character string containing the message to log.
#'
#' @return None
#' @export
log_message <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
}
