make_base_samples <- function() {
  tibble::tibble(
    sample = c("s1", "s2"),
    seqrun = "r1",
    fastq_R1 = c("a_R1.fastq.gz", "b_R1.fastq.gz"),
    fastq_R2 = c("a_R2.fastq.gz", "b_R2.fastq.gz"),
    neg_control = FALSE,
    pos_control = FALSE
  )
}

test_that("finalize_sample_table sets orient_group=orient for fwd-only", {
  withr::local_options(list(
    optimotu.pipeline.read_orientation = "fwd",
    optimotu.pipeline.do_rarefy = FALSE,
    optimotu.pipeline.duplicate_samples = "error"
  ))

  out <- optimotu.pipeline:::finalize_sample_table(make_base_samples())
  expect_true(all(out$orient == "fwd"))
  expect_true(all(out$orient_group == "fwd"))
  expect_equal(nrow(out), 2L)
})

test_that("finalize_sample_table marks all keys mix* when orient:mixed", {
  withr::local_options(list(
    optimotu.pipeline.read_orientation = "mixed",
    optimotu.pipeline.do_rarefy = FALSE,
    optimotu.pipeline.duplicate_samples = "error"
  ))

  out <- optimotu.pipeline:::finalize_sample_table(make_base_samples())
  expect_equal(sort(unique(out$orient)), c("fwd", "rev"))
  expect_equal(sort(unique(out$orient_group)), c("mixfwd", "mixrev"))
  expect_equal(nrow(out), 4L)
  expect_true(all(out$orient_group[out$orient == "fwd"] == "mixfwd"))
  expect_true(all(out$orient_group[out$orient == "rev"] == "mixrev"))
})

test_that("finalize_sample_table keeps disjoint custom orients as fwd/rev groups", {
  withr::local_options(list(
    optimotu.pipeline.read_orientation = "custom",
    optimotu.pipeline.do_rarefy = FALSE,
    optimotu.pipeline.duplicate_samples = "error",
    optimotu.pipeline.custom_sample_table = "dummy.tsv"
  ))

  samples <- make_base_samples() |>
    dplyr::mutate(orient = c("fwd", "rev"))
  out <- optimotu.pipeline:::finalize_sample_table(samples)
  expect_equal(out$orient_group, out$orient)
  expect_equal(nrow(out), 2L)
})

test_that("finalize_sample_table marks shared sample_keys as mixfwd/mixrev", {
  withr::local_options(list(
    optimotu.pipeline.read_orientation = "custom",
    optimotu.pipeline.do_rarefy = FALSE,
    optimotu.pipeline.duplicate_samples = "error",
    optimotu.pipeline.custom_sample_table = "dummy.tsv"
  ))

  samples <- tibble::tibble(
    sample = c("mixed", "mixed", "fwd_only"),
    seqrun = "r1",
    fastq_R1 = paste0(c("m1", "m2", "f1"), "_R1.fastq.gz"),
    fastq_R2 = paste0(c("m1", "m2", "f1"), "_R2.fastq.gz"),
    neg_control = FALSE,
    pos_control = FALSE,
    orient = c("fwd", "rev", "fwd")
  )
  out <- optimotu.pipeline:::finalize_sample_table(samples)
  mixed <- out[out$sample == "mixed", ]
  expect_equal(sort(mixed$orient_group), c("mixfwd", "mixrev"))
  expect_equal(out$orient_group[out$sample == "fwd_only"], "fwd")
})
