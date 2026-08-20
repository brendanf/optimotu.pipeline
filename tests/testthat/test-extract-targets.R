test_that("extract_targets uses runtime meta without calling tar_meta()", {
  rt <- targets::tar_runtime_object()
  old_meta <- rt$meta
  old_store <- rt$store
  old_target <- rt$target
  on.exit(
    {
      rt$meta <- old_meta
      rt$store <- old_store
      rt$target <- old_target
    },
    add = TRUE
  )
  rt$store <- NULL
  rt$target <- NULL
  rt$meta <- list(
    exists_record = function(x) x %in% c("stem1", "pat"),
    get_record = function(x) {
      switch(
        x,
        stem1 = list(type = "stem", name = "stem1"),
        pat = list(type = "pattern", children = c("pat_1", "pat_2"))
      )
    }
  )
  expect_equal(
    sort(optimotu.pipeline:::extract_targets(quote(list(stem1, pat, other)))),
    sort(c("stem1", "pat_1", "pat_2"))
  )
})

test_that("extract_targets uses the worker subpipeline when meta is unset", {
  rt <- targets::tar_runtime_object()
  old_meta <- rt$meta
  old_store <- rt$store
  old_target <- rt$target
  on.exit(
    {
      rt$meta <- old_meta
      rt$store <- old_store
      rt$target <- old_target
    },
    add = TRUE
  )
  rt$meta <- NULL
  rt$store <- tempfile("fake_store_")
  targets_env <- new.env(parent = emptyenv())
  pipeline <- list(targets = targets_env)
  stem <- new.env(parent = emptyenv())
  class(stem) <- c("tar_stem", "tar_builder", "tar_target")
  stem$name <- "stem1"
  targets_env$stem1 <- stem
  pat <- new.env(parent = emptyenv())
  class(pat) <- c("tar_pattern", "tar_target")
  pat$name <- "pat"
  pat$junction <- list(index = c(pat_1 = 1L, pat_2 = 2L))
  targets_env$pat <- pat
  target <- new.env(parent = emptyenv())
  target$subpipeline <- pipeline
  rt$target <- target
  expect_equal(
    sort(optimotu.pipeline:::extract_targets(quote(list(stem1, pat, other)))),
    sort(c("stem1", "pat_1", "pat_2"))
  )
})

test_that("extract_targets errors in-pipeline instead of calling tar_meta()", {
  rt <- targets::tar_runtime_object()
  old_meta <- rt$meta
  old_store <- rt$store
  old_target <- rt$target
  on.exit(
    {
      rt$meta <- old_meta
      rt$store <- old_store
      rt$target <- old_target
    },
    add = TRUE
  )
  rt$meta <- NULL
  rt$store <- tempfile("fake_store_")
  rt$target <- new.env(parent = emptyenv())
  expect_error(
    optimotu.pipeline:::extract_targets(quote(stem1)),
    "Cannot resolve target names"
  )
})

test_that("extract_targets does not call tar_meta() for an empty subpipeline", {
  rt <- targets::tar_runtime_object()
  old_meta <- rt$meta
  old_store <- rt$store
  old_target <- rt$target
  on.exit(
    {
      rt$meta <- old_meta
      rt$store <- old_store
      rt$target <- old_target
    },
    add = TRUE
  )
  rt$meta <- NULL
  rt$store <- "_targets"
  target <- new.env(parent = emptyenv())
  target$subpipeline <- list(targets = new.env(parent = emptyenv()))
  rt$target <- target
  expect_equal(
    length(optimotu.pipeline:::extract_targets(quote(stem1))),
    0L
  )
})
