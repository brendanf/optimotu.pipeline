#### Utility functions ####

#' Determine if the R process seems to be running inside a SLURM job
#'
#' @return a logical value indicating whether the R process is running inside a
#' SLURM job
is_slurm <- function() {
  nchar(Sys.getenv("SLURM_JOB_ID")) > 0 || nchar(Sys.which("sbatch")) > 0
}

#' Determine if the R process seems to be running "locally" (i.e. not on a cluster)
#' @return a logical value indicating whether the R process is running locally
is_local <- function() !is_slurm()

#' Determine if the R process seems to be running inside a Snakemake job
#' @return a logical value indicating whether the R process is running inside a
#' Snakemake job
is_snakemake <- function() !interactive() && exists("snakemake")

#' Get the number of CPUs available to the R process
#'
#' If running inside a Snakemake or SLURM job, this function will return the
#' number of CPUs allocated to the job. Otherwise, it will return the value of
#' the option `"optimotu_num_threads"` if it is set, otherwise the number of
#' CPUs available to the R process minus one.
#'
#' @return an integer indicating the number of CPUs available to the R process
#' @export
local_cpus <- function() {
  if (is_snakemake()) {
    snakemake@threads
  } else if (is_slurm()) {
    out <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
    # if slurm doesn't know how many cores, that means we're probably on the
    # login node, so we should only use 1 core.
    if (!checkmate::test_count(out, positive = TRUE)) out <- 1L
    out
  } else {
    getOption("optimotu_num_threads", max(parallel::detectCores() - 1L, 1L))
  }
}



#' Ensure that a table has at least one row
#'
#' @param table (`data.frame`) a table to check
#' @return (`data.frame`) the input table if it has at least one row, otherwise
#' a single-row table with the same columns and types, but all values `NA`. The
#' exception is if the table includes a column named "tar_group" (as produced by
#' [targets::tar_group()]), in which case that column is given the value `1L` in
#' the output if the input had 0 rows.
#' @keywords internal
#' @export
ensure_table <- function(table) {
  checkmate::assert_data_frame(table)
  if (nrow(table) > 0) {
    return(table)
  }
  out <- tibble::tibble(.rows = 1)
  for (n in names(table)) {
    if (n == "tar_group") {
      out[[n]] <- 1L
    } else if (is.factor(table[[n]])) {
      out[[n]] <- factor(NA, levels = levels(table[[n]]))
    } else if (is.integer(table[[n]])) {
      out[[n]] <- NA_integer_
    } else if (is.numeric(table[[n]])) {
      out[[n]] <- NA_real_
    } else if (is.character(table[[n]])) {
      out[[n]] <- NA_character_
    } else if (is.factor(table[[n]])) {
      out[[n]] <- factor(NA, levels = levels(table[[n]]))
    } else if (is.logical(table[[n]])) {
      out[[n]] <- NA
    } else if (is.list(table[[n]])) {
      out[[n]] <- list(NA)
    } else {
      stop("Unsupported column type: ", class(table[[n]]))
    }
  }
  out
}
