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
