#' Global ESS Result Methods
#'
#' @param x A \code{geoess_global_result} object.
#' @param ... Unused.
#'
#' @return The object, invisibly for \code{print}.
#' @export
print.geoess_global_result <- function(x, ...) {
  cat("geoess_global_result\n")
  cat("  type:", x$type, "\n")
  cat("  method:", x$method, "\n")
  cat("  n:", x$n, "\n")
  cat("  estimate:\n")
  print(x$estimate)
  if (!is.null(x$mcse)) {
    cat("  mcse:\n")
    print(x$mcse)
  }
  if (!is.null(x$weights_summary)) {
    cat("  weights:\n")
    print(x$weights_summary)
  }
  invisible(x)
}

#' @export
summary.geoess_global_result <- function(object, ...) {
  object
}
