#' Distribute tasks to bins
#' 
#' This function distributes tasks to bins using a greedy algorithm.
#' Each task is assigned to the first bin that has enough space;
#' or if no such bin exists, a new bin is created. Tasks with weight
#' greater than the maximum bin weight are assigned to a bin of their
#' own.
#'
#' @param w (`numeric`) weight of each task
#' @param max_w (`numeric` scalar) maximum bin weight
#'
#' @return (`integer` vector) bin assignments
#' @export
#'
#' @examples
#' w <- rev(1:8)
#' max_w <- 5
#' distribute_tasks(w, max_w)
#' # [1] 1 2 3 4 5 6 6 5

distribute_tasks <- function(w, max_w) {
  checkmate::assert_numeric(w, lower = 0, finite = TRUE)
  checkmate::assert_number(max_w, lower = 0, finite = TRUE)

  # weight assigned to each bin
  bins <- numeric()
  # bin assignment for each task
  out <- integer(length(w))
  for (i in seq_along(w)) {
    if (w[i] > max_w) {
      bins <- c(bins, w[i])
      out[i] <- length(bins)
    } else {
      found <- FALSE
      for (j in seq_along(bins)) {
        if (w[i] + bins[j] <= max_w) {
          out[i] <- j
          bins[j] <- bins[j] + w[i]
          found <- TRUE
          break
        }
      }
      if (!found) {
        bins <- c(bins, w[i])
        out[i] <- length(bins)
      }
    }
  }
  out
}
