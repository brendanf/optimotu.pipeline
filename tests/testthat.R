# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(optimotu.pipeline)

# lulu::lulu() (called from test_lulu.R) unconditionally writes a verbose
# "lulu.log_<timestamp>" file to the working directory, with no option to
# suppress or redirect it. These logs are useful for debugging failures in
# our own re-implementation, so we only delete the ones created by this run
# if the whole suite passes; if anything fails, they are left behind
# alongside the (still-propagated) test failure for post-mortem inspection.
old_logs <- list.files(pattern = "^lulu\\.log_")
tests_passed <- TRUE

tryCatch(
  test_check("optimotu.pipeline"),
  error = function(e) {
    tests_passed <<- FALSE
    stop(e)
  },
  finally = {
    new_logs <- setdiff(list.files(pattern = "^lulu\\.log_"), old_logs)
    if (tests_passed) {
      unlink(new_logs)
    } else if (length(new_logs)) {
      message(
        "Tests failed; keeping LULU debug log(s) for inspection: ",
        paste(new_logs, collapse = ", ")
      )
    }
  }
)
