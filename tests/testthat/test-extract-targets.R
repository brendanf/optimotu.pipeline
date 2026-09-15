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
  # junction_init() builds the version-appropriate junction shape
  # ($index names in targets <= 1.10; $splits in >= 1.11).
  pat$junction <- targets:::junction_init(
    nexus = "pat",
    splits = c("pat_1", "pat_2")
  )
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

test_that("extract_targets resolves dynamic co-branch builders by stem", {
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
  # Worker subpipeline for a map() branch often has only the sibling child,
  # not the parent pattern object named "lulu_match_table".
  child <- new.env(parent = emptyenv())
  class(child) <- c("tar_stem", "tar_builder", "tar_target")
  child$name <- "lulu_match_table_abc123def4567890"
  targets_env[[child$name]] <- child
  target <- new.env(parent = emptyenv())
  target$subpipeline <- pipeline
  rt$target <- target
  expect_equal(
    optimotu.pipeline:::extract_targets(quote(lulu_match_table)),
    "lulu_match_table_abc123def4567890"
  )
})

test_that("tar_stem strips dynamic branch hashes", {
  expect_equal(
    optimotu.pipeline:::tar_stem(
      c(
        "seqtable_raw_LIFEPLAN.00001_eebce0ec729d5008",
        "seqtable_raw_LIFEPLAN.00001"
      )
    ),
    c("seqtable_raw_LIFEPLAN.00001", "seqtable_raw_LIFEPLAN.00001")
  )
})
