# Skip helpers for optional external executables used in integration tests.

skip_if_no_vsearch <- function() {
  tc_skip_if_no_exec("vsearch")
}

skip_if_no_usearch <- function() {
  tc_skip_if_no_exec("usearch")
}
