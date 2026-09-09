test_that("unnest_yaml_list flattens length-1 nested lists", {
  nested <- list(list(a = 1), list(b = 2))
  out <- optimotu.pipeline:::unnest_yaml_list(nested)
  expect_equal(out, list(a = 1, b = 2))
})

test_that("parse_project_name handles default, valid, and invalid names", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  expect_warning(
    out <- optimotu.pipeline:::parse_project_name(list()),
    "Missing project name"
  )
  expect_equal(out, "metabarcoding_project")
  expect_equal(optimotu.pipeline::project_name(), "metabarcoding_project")

  optimotu.pipeline:::parse_project_name(list(project_name = "my_project-1"))
  expect_equal(optimotu.pipeline::project_name(), "my_project-1")

  expect_error(
    optimotu.pipeline:::parse_project_name(list(project_name = "bad name")),
    "Project name should consist"
  )
})

test_that("parse_orient sets default and explicit orientation", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  expect_message(
    optimotu.pipeline:::parse_orient(list(orient = NULL)),
    "No read orientation specified"
  )
  expect_equal(optimotu.pipeline::read_orientation(), "fwd")

  optimotu.pipeline:::parse_orient(list(orient = "mixed"))
  expect_equal(optimotu.pipeline::read_orientation(), "mixed")

  expect_error(
    optimotu.pipeline:::parse_orient(list(orient = "forward")),
    "Must comply to pattern"
  )
})

test_that("parse_custom_sample_table accepts FALSE and valid path", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  expect_null(optimotu.pipeline::custom_sample_table())
  optimotu.pipeline:::parse_custom_sample_table(
    list(custom_sample_table = FALSE)
  )
  expect_null(optimotu.pipeline::custom_sample_table())

  sample_tbl <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("sample\tseqrun", "s1\tr1"), sample_tbl)
  optimotu.pipeline:::parse_custom_sample_table(
    list(custom_sample_table = sample_tbl)
  )
  expect_equal(optimotu.pipeline::custom_sample_table(), sample_tbl)
  expect_true(optimotu.pipeline:::do_custom_sample_table())
})

test_that("parse_added_reference requires both files and sets options", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())
  options(
    optimotu.pipeline.do_protax = TRUE,
    optimotu.pipeline.do_added_reference = FALSE,
    optimotu.pipeline.added_reference_fasta = NULL,
    optimotu.pipeline.added_reference_table = NULL
  )

  fasta <- withr::local_tempfile(fileext = ".fa")
  table <- withr::local_tempfile(fileext = ".xlsx")
  writeLines(">x\nACGT", fasta)
  writeLines("placeholder", table)

  expect_error(
    suppressWarnings(
      optimotu.pipeline:::parse_added_reference(
        list(added_reference = list(fasta = fasta))
      )
    ),
    "both must be given"
  )

  expect_warning(
    optimotu.pipeline:::parse_added_reference(
      list(added_reference = list(list(fasta = fasta), list(table = table)))
    ),
    "deprecated"
  )
  expect_true(optimotu.pipeline::do_added_reference())
  expect_equal(optimotu.pipeline::added_reference_fasta(), fasta)
  expect_equal(optimotu.pipeline::added_reference_table(), table)
})

test_that("parse_supplemental_asv_options parses enabled sets", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  fasta1 <- withr::local_tempfile(fileext = ".fasta")
  fasta2 <- withr::local_tempfile(fileext = ".fasta")
  taxfile <- withr::local_tempfile(fileext = ".tsv")
  sample_tbl <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c(">suppA1", "ACGT"), fasta1)
  writeLines(c(">suppB1", "TGCA"), fasta2)
  writeLines("seq_id\tkingdom\nsuppB1\tFungi", taxfile)
  writeLines("sample\tseq_id\tnread\ns1\tsuppA1\t2", sample_tbl)

  optimotu.pipeline:::parse_supplemental_asv_options(
    list(
      supplemental_asv = list(
        list(
          env = list(
            sequences = fasta1,
            sample_table = sample_tbl,
            taxonomy = FALSE
          )
        ),
        list(
          ref = list(
            sequences = fasta2,
            taxonomy = taxfile,
            enabled = TRUE
          )
        )
      )
    )
  )

  expect_true(optimotu.pipeline::do_supp_asv())
  expect_setequal(
    optimotu.pipeline::supp_asv_set_names(),
    c("env", "ref")
  )
  expect_equal(optimotu.pipeline::supp_asv_sequences("env"), fasta1)
  expect_equal(
    optimotu.pipeline::supp_asv_sample_table("env"),
    sample_tbl
  )
  expect_equal(optimotu.pipeline::supp_asv_taxonomy_mode("env"), "none")
  expect_equal(optimotu.pipeline::supp_asv_taxonomy_mode("ref"), "file")
  expect_equal(
    optimotu.pipeline::supp_asv_taxonomy_file("ref"),
    taxfile
  )
})

test_that("parse_supplemental_asv_options supports header taxonomy", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  fasta <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">id1;tax=d:Eukaryota,k:Fungi,s:sp", "ACGT"), fasta)

  optimotu.pipeline:::parse_supplemental_asv_options(
    list(
      supplemental_asv = list(
        hdr = list(
          sequences = fasta,
          taxonomy = TRUE
        )
      )
    )
  )

  expect_true(optimotu.pipeline::do_supp_asv())
  expect_equal(
    optimotu.pipeline::supp_asv_taxonomy_mode("hdr"),
    "header"
  )
  expect_null(optimotu.pipeline::supp_asv_taxonomy_file("hdr"))
})

test_that("parse_supplemental_asv_options rejects whitespace in headers", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  fasta <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">bad id", "ACGT"), fasta)

  expect_error(
    optimotu.pipeline:::parse_supplemental_asv_options(
      list(
        supplemental_asv = list(
          bad = list(sequences = fasta)
        )
      )
    ),
    "must not contain whitespace or semicolons"
  )
})

test_that("parse_supplemental_asv_options rejects semicolons in headers", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  fasta <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">bad;id", "ACGT"), fasta)

  expect_error(
    optimotu.pipeline:::parse_supplemental_asv_options(
      list(
        supplemental_asv = list(
          bad = list(sequences = fasta)
        )
      )
    ),
    "must not contain whitespace or semicolons"
  )
})

test_that("parse_supplemental_asv_options allows semicolons with header taxonomy", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  fasta <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">id1;tax=d:Eukaryota,k:Fungi,s:sp", "ACGT"), fasta)

  optimotu.pipeline:::parse_supplemental_asv_options(
    list(
      supplemental_asv = list(
        hdr = list(
          sequences = fasta,
          taxonomy = TRUE
        )
      )
    )
  )

  expect_equal(
    optimotu.pipeline::supp_asv_taxonomy_mode("hdr"),
    "header"
  )
})

test_that("supplemental ASV accessors accept vector set names", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  fasta1 <- withr::local_tempfile(fileext = ".fasta")
  fasta2 <- withr::local_tempfile(fileext = ".fasta")
  sample_tbl <- withr::local_tempfile(fileext = ".tsv")
  taxfile <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c(">env1", "ACGT"), fasta1)
  writeLines(c(">ref1", "TGCA"), fasta2)
  writeLines(c("sample\tseq_id\tnread", "s1\tenv1\t5"), sample_tbl)
  writeLines(
    c(
      "seq_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies",
      "ref1\tEukaryota\tChordata\tActinopteri\tPerciformes\tSerranidae\tEpinephelus\tcoioides"
    ),
    taxfile
  )

  optimotu.pipeline:::parse_supplemental_asv_options(
    list(
      supplemental_asv = list(
        env = list(sequences = fasta1, sample_table = sample_tbl),
        ref = list(sequences = fasta2, taxonomy = taxfile)
      )
    )
  )

  expect_equal(
    optimotu.pipeline::supp_asv_sequences(c("ref", "env")),
    c(ref = fasta2, env = fasta1)
  )
  expect_equal(
    optimotu.pipeline::supp_asv_sample_table(c("env", "ref")),
    c(env = sample_tbl, ref = NA_character_)
  )
  expect_equal(
    optimotu.pipeline::supp_asv_taxonomy_mode(c("env", "ref")),
    c(env = "none", ref = "file")
  )
  expect_equal(
    optimotu.pipeline::supp_asv_taxonomy_file(c("env", "ref")),
    c(env = NA_character_, ref = taxfile)
  )
  expect_equal(
    names(optimotu.pipeline::supp_asv_set(c("ref", "env"))),
    c("ref", "env")
  )
  expect_error(
    optimotu.pipeline::supp_asv_set(c("env", "missing")),
    "Unknown supplemental ASV set\\(s\\)"
  )
})

test_that("parse_supplemental_asv_options rejects duplicate set names early", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  expect_error(
    optimotu.pipeline:::parse_supplemental_asv_options(
      list(
        supplemental_asv = list(
          list(dup = list(sequences = "does_not_exist_1.fasta")),
          list(dup = list(sequences = "does_not_exist_2.fasta"))
        )
      )
    ),
    "Duplicate supplemental ASV set name\\(s\\)"
  )
})

test_that("parse_supplemental_asv_options allows duplicate source IDs", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  fasta1 <- withr::local_tempfile(fileext = ".fasta")
  fasta2 <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">dup1", "ACGT"), fasta1)
  writeLines(c(">dup1", "TGCA"), fasta2)

  optimotu.pipeline:::parse_supplemental_asv_options(
    list(
      supplemental_asv = list(
        a = list(sequences = fasta1),
        b = list(sequences = fasta2)
      )
    )
  )
  expect_setequal(optimotu.pipeline::supp_asv_set_names(), c("a", "b"))
})

test_that("parse_parallel_options handles workers precedence and bounds", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  expect_warning(
    optimotu.pipeline:::parse_parallel_options(
      list(
        local_threads = 3L,
        max_batchsize = 1000L,
        workers_per_seqrun = 4L,
        jobs_per_seqrun = 2L,
        min_workers = 2L,
        max_workers = 10L
      )
    ),
    "both 'workers_per_seqrun' and 'jobs_per_seqrun'"
  )

  expect_equal(getOption("optimotu_num_threads"), 3L)
  expect_equal(optimotu.pipeline::max_batchsize(), 1000L)
  expect_equal(optimotu.pipeline::workers_per_seqrun(), 4L)
  expect_equal(optimotu.pipeline::min_workers(), 2L)
  expect_equal(optimotu.pipeline::max_workers(), 10L)
})

test_that("parse_rarefy_options supports number and fraction modes", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  optimotu.pipeline:::parse_rarefy_options(
    list(rarefy = list(number = 100L))
  )
  expect_equal(optimotu.pipeline::rarefy_number(), 100L)
  expect_true(optimotu.pipeline::do_rarefy())

  optimotu.pipeline:::parse_rarefy_options(
    list(rarefy = list(numerator = 1L, denominator = 10L))
  )
  expect_equal(optimotu.pipeline::rarefy_numerator(), 1L)
  expect_equal(optimotu.pipeline::rarefy_denominator(), 10L)

  expect_error(
    optimotu.pipeline:::parse_rarefy_options(
      list(rarefy = list(number = 100L, numerator = 1L))
    ),
    "cannot be used in conjunction"
  )
})

test_that("parse_denoising_options defaults to dada2 and accepts unoise", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  optimotu.pipeline:::parse_denoising_options(list())
  expect_equal(optimotu.pipeline::denoising_method(), "dada2")
  expect_true(optimotu.pipeline::do_dada2())
  expect_false(optimotu.pipeline::do_unoise())
  expect_equal(optimotu.pipeline::denoising_pool(), "sample")

  optimotu.pipeline:::parse_denoising_options(list(denoising = "unoise"))
  expect_equal(optimotu.pipeline::denoising_method(), "unoise")
  expect_true(optimotu.pipeline::do_unoise())
  expect_equal(optimotu.pipeline::unoise_alpha(), 2)
  expect_equal(optimotu.pipeline::unoise_minsize(), 8L)

  optimotu.pipeline:::parse_denoising_options(
    list(
      denoising = list(
        method = "unoise",
        pool = "sample",
        unoise = list(
          alpha = 1.5,
          minsize = 4L
        )
      )
    )
  )
  expect_equal(optimotu.pipeline::unoise_alpha(), 1.5)
  expect_equal(optimotu.pipeline::unoise_minsize(), 4L)

  expect_error(
    optimotu.pipeline:::parse_denoising_options(list(denoising = "nope")),
    "Must be element of set"
  )
  expect_error(
    optimotu.pipeline:::parse_denoising_options(
      list(denoising = list(method = "poanoise"))
    ),
    "not implemented"
  )
  expect_error(
    optimotu.pipeline:::parse_denoising_options(
      list(denoising = list(method = "unoise", pool = "project"))
    ),
    "only 'sample' is implemented"
  )
  expect_error(
    optimotu.pipeline:::parse_denoising_options(
      list(
        denoising = list(
          method = "unoise",
          unoise = list(merge = list(min_overlap = 20L))
        )
      )
    ),
    "moved to the top-level 'merging:'"
  )
})

test_that("parse_merge_options uses method-specific defaults", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  options(
    optimotu.pipeline.denoising_method = "dada2",
    optimotu.pipeline.merge_min_overlap = NULL,
    optimotu.pipeline.merge_max_mismatch = NULL
  )
  expect_equal(optimotu.pipeline::merge_min_overlap(), 10L)
  expect_equal(optimotu.pipeline::merge_max_mismatch(), 1)

  options(optimotu.pipeline.denoising_method = "unoise")
  expect_equal(optimotu.pipeline::merge_min_overlap(), 16L)
  expect_equal(optimotu.pipeline::merge_max_mismatch(), 5)

  options(optimotu.pipeline.denoising_method = "dada2")
  optimotu.pipeline:::parse_merge_options(
    list(merging = list(min_overlap = 12L, max_mismatch = 2))
  )
  expect_equal(optimotu.pipeline::merge_min_overlap(), 12L)
  expect_equal(optimotu.pipeline::merge_max_mismatch(), 2)

  expect_error(
    optimotu.pipeline:::parse_merge_options(
      list(merging = list(max_mismatch = 0.1))
    ),
    "integer count"
  )

  options(optimotu.pipeline.denoising_method = "unoise")
  options(optimotu.pipeline.merge_max_mismatch = NULL)
  optimotu.pipeline:::parse_merge_options(
    list(merging = list(max_mismatch = 0.05))
  )
  expect_equal(optimotu.pipeline::merge_max_mismatch(), 0.05)
})

test_that("merge_dist_config_lists inherits same-method keys only", {
  expect_equal(
    optimotu.pipeline:::merge_dist_config_lists(NULL, NULL),
    list(method = "usearch")
  )
  expect_equal(
    optimotu.pipeline:::merge_dist_config_lists(
      NULL,
      list(method = "hamming")
    ),
    list(method = "hamming")
  )
  expect_equal(
    optimotu.pipeline:::merge_dist_config_lists(
      list(method = "usearch", usearch = "bin/usearch"),
      list(method = "usearch")
    ),
    list(method = "usearch", usearch = "bin/usearch")
  )
  expect_equal(
    optimotu.pipeline:::merge_dist_config_lists(
      list(method = "wfa2", match = 0),
      list(method = "usearch", usearch = "bin/usearch")
    ),
    list(method = "wfa2", match = 0)
  )
  expect_equal(
    optimotu.pipeline:::merge_dist_config_lists(
      list(usearch = "bin/usearch"),
      list(method = "usearch")
    ),
    list(method = "usearch", usearch = "bin/usearch")
  )
})

test_that("parse_executable_options stores configured paths", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())
  options(optimotu.pipeline.executables = NULL)

  optimotu.pipeline:::parse_executable_options(
    list(executables = list(usearch = "bin/usearch", vsearch = "vsearch"))
  )
  expect_equal(
    optimotu.pipeline::configured_executables(),
    c(usearch = "bin/usearch", vsearch = "vsearch")
  )
})

test_that("added_reference under taxonomy.protax and deprecated top-level", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())
  options(
    optimotu.pipeline.do_added_reference = FALSE,
    optimotu.pipeline.added_reference_fasta = NULL,
    optimotu.pipeline.added_reference_table = NULL,
    optimotu.pipeline.do_protax = FALSE
  )

  fasta <- withr::local_tempfile(fileext = ".fasta")
  table <- withr::local_tempfile(fileext = ".xlsx")
  writeLines(">s1\nACGT", fasta)
  # empty file is enough for assert_file_exists
  file.create(table)

  optimotu.pipeline:::parse_protax_options(
    list(
      aligned = FALSE,
      location = tempdir(),
      added_reference = list(fasta = fasta, table = table)
    )
  )
  expect_true(optimotu.pipeline::do_added_reference())
  expect_equal(optimotu.pipeline::added_reference_fasta(), fasta)
  expect_equal(optimotu.pipeline::added_reference_table(), table)

  options(
    optimotu.pipeline.do_added_reference = FALSE,
    optimotu.pipeline.added_reference_fasta = NULL,
    optimotu.pipeline.added_reference_table = NULL,
    optimotu.pipeline.do_protax = TRUE
  )
  expect_warning(
    optimotu.pipeline:::parse_added_reference(
      list(added_reference = list(fasta = fasta, table = table))
    ),
    "deprecated"
  )
  expect_true(optimotu.pipeline::do_added_reference())

  options(
    optimotu.pipeline.do_added_reference = FALSE,
    optimotu.pipeline.do_protax = FALSE
  )
  expect_error(
    optimotu.pipeline:::parse_added_reference(
      list(added_reference = list(fasta = fasta, table = table))
    ),
    "only valid with the Protax classifier"
  )

  # empty stub is ignored
  expect_silent(
    optimotu.pipeline:::parse_added_reference(
      list(added_reference = list(fasta = NULL, table = NULL))
    )
  )
})

test_that("parse_filter_options is method-aware for paired vs merged keys", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())
  options(optimotu.pipeline.merged_filter_options = NULL)

  options(optimotu.pipeline.denoising_method = "dada2")
  expect_warning(
    optimotu.pipeline:::parse_filter_options(
      list(filtering = list(maxEE_R1 = 1, maxEE_R2 = 2, maxEE = 0.5))
    ),
    "Ignoring merged-read filtering"
  )
  expect_equal(unname(optimotu.pipeline::dada2_maxEE()["maxEE_R1"]), 1)

  options(optimotu.pipeline.denoising_method = "unoise")
  expect_warning(
    optimotu.pipeline:::parse_filter_options(
      list(filtering = list(maxEE = 0.8, minLen = 100L, maxEE_R1 = 2))
    ),
    "Ignoring paired-read filtering"
  )
  opts <- merged_filter_options()
  expect_equal(opts$maxEE, 0.8)
  expect_equal(opts$minLen, 100L)
})

test_that("parse_cluster_options reads min_parallel_ops and max_batch_ops", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  expect_equal(
    optimotu.pipeline:::cluster_ops_defaults("hamming"),
    list(min_parallel_ops = 1e6, max_batch_ops = 1e10)
  )
  expect_equal(
    optimotu.pipeline:::cluster_ops_defaults("usearch"),
    list(min_parallel_ops = 1e6, max_batch_ops = 1e10)
  )
  expect_equal(
    optimotu.pipeline:::cluster_ops_defaults("wfa2"),
    list(min_parallel_ops = 1e4, max_batch_ops = 1e8)
  )
  expect_equal(
    optimotu.pipeline:::cluster_ops_defaults("edlib"),
    list(min_parallel_ops = 1e4, max_batch_ops = 1e8)
  )
  expect_equal(
    optimotu.pipeline:::cluster_ops_defaults("ksw2"),
    list(min_parallel_ops = 1e4, max_batch_ops = 1e8)
  )
  expect_equal(
    optimotu.pipeline:::cluster_ops_defaults("hybrid"),
    list(min_parallel_ops = 1e4, max_batch_ops = 1e8)
  )

  msgs <- testthat::capture_messages(
    optimotu.pipeline:::parse_cluster_options(
      list(clustering = list(dist_config = "wfa2"))
    )
  )
  expect_true(any(grepl("min_parallel_ops", msgs)))
  expect_true(any(grepl("max_batch_ops", msgs)))
  expect_equal(cluster_min_parallel_ops(), 1e4)
  expect_equal(cluster_max_batch_ops(), 1e8)

  msgs <- testthat::capture_messages(
    optimotu.pipeline:::parse_cluster_options(
      list(
        clustering = list(
          dist_config = "wfa2",
          min_parallel_ops = 2e6,
          max_batch_ops = 5e9
        )
      )
    )
  )
  expect_false(any(grepl("job sizing", msgs)))
  expect_equal(cluster_min_parallel_ops(), 2e6)
  expect_equal(cluster_max_batch_ops(), 5e9)

  expect_error(
    suppressMessages(
      optimotu.pipeline:::parse_cluster_options(
        list(
          clustering = list(
            dist_config = "wfa2",
            min_parallel_ops = 1e8,
            max_batch_ops = 1e6
          )
        )
      )
    ),
    "max_batch_ops must be greater than or equal to"
  )
})

test_that("parse_cluster_thresholds handles self, file, and load-only modes", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  train_fa <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">s1;tax=k:Fungi,p:Asco", "ACGT"), train_fa)

  optimotu.pipeline:::parse_cluster_thresholds(
    list(train_data = "self", min_conf = 0.6, dist_max = 0.3)
  )
  expect_true(do_optimize_thresholds())
  expect_true(do_optimize_thresholds_self())
  expect_false(do_optimize_thresholds_reference())
  expect_false(do_optimize_thresholds_file())
  expect_equal(cluster_min_conf(), 0.6)
  expect_equal(cluster_dist_max(), 0.3)

  optimotu.pipeline:::parse_cluster_thresholds(
    list(train_data = "reference", min_taxa = 3L)
  )
  expect_true(do_optimize_thresholds())
  expect_true(do_optimize_thresholds_reference())
  expect_false(do_optimize_thresholds_self())
  expect_false(do_optimize_thresholds_file())
  expect_equal(cluster_min_taxa(), 3)

  optimotu.pipeline:::parse_cluster_thresholds(
    list(train_data = train_fa, file = "output/optima.tsv")
  )
  expect_true(do_optimize_thresholds())
  expect_true(do_optimize_thresholds_file())
  expect_equal(optimize_thresholds_file(), train_fa)
  expect_equal(cluster_thresholds(), "output/optima.tsv")

  load_tsv <- withr::local_tempfile(fileext = ".tsv")
  writeLines("rank\tthreshold", load_tsv)
  options(
    optimotu.pipeline.do_optimize_thresholds = NULL,
    optimotu.pipeline.do_optimize_thresholds_self = NULL,
    optimotu.pipeline.do_optimize_thresholds_reference = NULL,
    optimotu.pipeline.do_optimize_thresholds_file = NULL,
    optimotu.pipeline.optimize_thresholds_file = NULL
  )
  optimotu.pipeline:::parse_cluster_thresholds(load_tsv)
  expect_false(isTRUE(do_optimize_thresholds()))
  expect_equal(cluster_thresholds(), load_tsv)

  options(
    optimotu.pipeline.do_optimize_thresholds = NULL,
    optimotu.pipeline.clustering_thresholds = NULL
  )
  optimotu.pipeline:::parse_cluster_thresholds(list(file = load_tsv))
  expect_false(isTRUE(do_optimize_thresholds()))
  expect_equal(cluster_thresholds(), load_tsv)
})

test_that("parse_cluster_options resolves memory_budget_mb defaults", {
  old <- options()
  withr::defer(options(old), testthat::teardown_env())

  suppressMessages(
    optimotu.pipeline:::parse_cluster_options(
      list(
        clustering = list(
          thresholds = list(train_data = "self"),
          dist_config = "wfa2"
        )
      )
    )
  )
  expect_identical(cluster_memory_budget_mb(), "auto")

  expect_identical(
    optimotu.pipeline:::resolve_cluster_memory_budget_mb(
      NULL,
      method = "ksw2",
      do_optimize = TRUE
    ),
    "auto"
  )

  suppressMessages(
    optimotu.pipeline:::parse_cluster_options(
      list(
        clustering = list(
          thresholds = list(train_data = "self"),
          dist_config = "usearch"
        )
      )
    )
  )
  expect_null(cluster_memory_budget_mb())

  suppressMessages(
    optimotu.pipeline:::parse_cluster_options(
      list(
        clustering = list(
          thresholds = list(train_data = "self"),
          dist_config = "wfa2",
          memory_budget_mb = 256
        )
      )
    )
  )
  expect_equal(cluster_memory_budget_mb(), 256)

  expect_error(
    suppressMessages(
      optimotu.pipeline:::parse_cluster_options(
        list(
          clustering = list(
            thresholds = list(train_data = "self"),
            dist_config = "usearch",
            memory_budget_mb = "auto"
          )
        )
      )
    ),
    "not supported with dist_config method 'usearch'"
  )

  # Load-only thresholds: no auto budget even for native methods
  load_tsv <- withr::local_tempfile(fileext = ".tsv")
  writeLines("rank\tthreshold", load_tsv)
  suppressMessages(
    optimotu.pipeline:::parse_cluster_options(
      list(
        clustering = list(
          thresholds = load_tsv,
          dist_config = "wfa2"
        )
      )
    )
  )
  expect_null(cluster_memory_budget_mb())
})

test_that("cluster_clust_config and cluster_parallel_config return calls", {
  expect_identical(
    cluster_clust_config("hamming"),
    quote(optimotu::clust_slink())
  )
  expect_identical(
    cluster_clust_config("wfa2"),
    quote(optimotu::clust_tree())
  )
  expect_identical(
    cluster_clust_config("ksw2"),
    quote(optimotu::clust_tree())
  )
  expect_identical(
    cluster_parallel_config(4L, "hamming"),
    quote(optimotu::parallel_merge(threads = 4L))
  )
  expect_identical(
    cluster_parallel_config(2L, "usearch"),
    quote(optimotu::parallel_concurrent(threads = 2L))
  )
  expect_identical(
    cluster_parallel_config(2L, "ksw2"),
    quote(optimotu::parallel_concurrent(threads = 2L))
  )
  expect_identical(
    cluster_parallel_config(local_cpus(), "hamming"),
    quote(optimotu::parallel_merge(threads = local_cpus()))
  )
})
