test_that("synthetic taxonomic fixture has expected hierarchy and split", {
  x <- tc_generate_dataset(seed = 123L)
  tc_assert_dataset_shape(x)
})

test_that("sintax runs on generated reference/query fixtures", {
  tc_skip_if_no_pkg("Biostrings")
  vsearch <- tc_skip_if_no_exec("vsearch")
  x <- tc_generate_dataset(seed = 123L)
  tc_assert_dataset_shape(x)

  ref_file <- withr::local_tempfile(fileext = ".fasta")
  query_file <- withr::local_tempfile(fileext = ".fasta")
  tc_write_sintax_reference(x$reference, ref_file)
  tc_write_plain_query(x$query, query_file)

  old <- options(VSEARCH = vsearch, vsearch = vsearch)
  withr::defer(options(old), testthat::teardown_env())

  out <- optimotu.pipeline::sintax(
    query = query_file,
    ref = ref_file,
    ncpu = 1L,
    id_is_int = TRUE
  )

  checkmate::expect_names(
    names(out),
    must.include = c("seq_idx", "rank", "parent_taxonomy", "taxon", "prob")
  )
  testthat::expect_true(is.integer(out$seq_idx))
  p <- as.numeric(out$prob)
  testthat::expect_true(any(!is.na(p)))
  testthat::expect_true(all(p[!is.na(p)] >= 0))
  testthat::expect_true(all(p[!is.na(p)] <= 1))
})

test_that("epa-ng and gappa run on generated fixtures", {
  tc_skip_if_no_pkg("Biostrings")
  epa_exec <- tc_skip_if_no_exec("epa-ng")
  gappa_exec <- {
    p <- tc_find_gappa()
    if (!nzchar(p)) {
      testthat::skip("Missing executable: gappa")
    }
    p
  }
  iqtree_exec <- tc_skip_if_no_exec(c("iqtree3", "iqtree"), "iqtree3 or iqtree")

  x <- tc_generate_dataset(seed = 123L)
  tc_assert_dataset_shape(x)

  ref_file <- withr::local_tempfile(fileext = ".fasta")
  query_file <- withr::local_tempfile(fileext = ".fasta")
  tc_write_plain_query(x$reference, ref_file)
  tc_write_plain_query(x$query, query_file)

  prefix <- file.path(withr::local_tempdir(), "tiny_iqtree")
  iq <- processx::run(
    iqtree_exec,
    c("-s", ref_file, "-m", "GTR+F", "-nt", "1", "-redo", "-pre", prefix),
    stderr = "|",
    stdout = "|",
    error_on_status = FALSE
  )
  testthat::expect_equal(iq$status, 0L)

  tree_file <- paste0(prefix, ".treefile")
  model_file <- paste0(prefix, ".iqtree")
  testthat::expect_true(file.exists(tree_file))
  testthat::expect_true(file.exists(model_file))

  model <- optimotu.pipeline:::parse_iqtree_model(file = model_file)
  testthat::expect_true(nchar(model) > 0)

  old_gappa <- Sys.getenv("gappa", unset = NA_character_)
  Sys.setenv(gappa = gappa_exec)
  withr::defer(
    {
      if (is.na(old_gappa)) {
        Sys.unsetenv("gappa")
      } else {
        Sys.setenv(gappa = old_gappa)
      }
    },
    testthat::teardown_env()
  )

  old <- options(
    `epa-ng` = epa_exec,
    EPA_NG = epa_exec,
    gappa = gappa_exec,
    GAPPA = gappa_exec
  )
  withr::defer(options(old), testthat::teardown_env())

  outdir <- withr::local_tempdir()
  jplace <- optimotu.pipeline::epa_ng(
    ref_msa = ref_file,
    tree = tree_file,
    query = query_file,
    outdir = outdir,
    model = model,
    ncpu = 1L,
    exec = epa_exec
  )

  testthat::expect_true(file.exists(jplace))
  jp <- jsonlite::read_json(jplace)
  testthat::expect_true(optimotu.pipeline:::is_jplace(jp))

  taxonomy_file <- withr::local_tempfile(fileext = ".tsv")
  tc_write_gappa_taxonomy(x$reference, taxonomy_file)

  out <- optimotu.pipeline::gappa_assign(
    jplace = jplace,
    taxonomy = taxonomy_file,
    outgroup = x$reference$seq_id[[1]],
    ranks = optimotu.pipeline::tax_ranks(),
    ncpu = 1L,
    allow_file_overwriting = TRUE,
    id_is_int = TRUE
  )

  checkmate::expect_names(
    names(out),
    must.include = c("seq_idx", "rank", "parent_taxon", "taxon", "prob")
  )
  testthat::expect_true(is.integer(out$seq_idx))
  testthat::expect_true(all(out$prob >= 0 & out$prob <= 1))
})

test_that("epa-ng direct IQ-TREE model-file parsing works across model matrix", {
  tc_skip_if_no_pkg("Biostrings")
  epa_exec <- tc_skip_if_no_exec("epa-ng")
  iqtree_exec <- tc_skip_if_no_exec(c("iqtree3", "iqtree"), "iqtree3 or iqtree")

  x <- tc_generate_dataset(seed = 123L)
  tc_assert_dataset_shape(x)

  ref_file <- withr::local_tempfile(fileext = ".fasta")
  query_file <- withr::local_tempfile(fileext = ".fasta")
  tc_write_plain_query(x$reference, ref_file)
  tc_write_plain_query(x$query, query_file)

  for (model_name in tc_iqtree_model_matrix()) {
    prefix <- file.path(
      withr::local_tempdir(),
      paste0("tiny_iqtree_", gsub("[^A-Za-z0-9]+", "_", model_name))
    )
    iq <- processx::run(
      iqtree_exec,
      c("-s", ref_file, "-m", model_name, "-nt", "1", "-redo", "-pre", prefix),
      stderr = "|",
      stdout = "|",
      error_on_status = FALSE
    )
    testthat::expect_equal(
      iq$status,
      0L,
      info = paste("IQ-TREE model", model_name)
    )

    tree_file <- paste0(prefix, ".treefile")
    model_file <- paste0(prefix, ".iqtree")
    testthat::expect_true(file.exists(tree_file), info = model_name)
    testthat::expect_true(file.exists(model_file), info = model_name)

    expected <- tc_parse_iqtree_expected(model_file)
    outdir <- withr::local_tempdir()
    epa <- tc_run_epa_capture(
      epa_exec = epa_exec,
      ref_file = ref_file,
      query_file = query_file,
      tree_file = tree_file,
      model_arg = model_file,
      outdir = outdir
    )
    epa_log <- paste(c(epa$stdout, epa$stderr), collapse = "\n")
    direct_ok <- epa$status == 0L

    parsed_model <- optimotu.pipeline:::parse_iqtree_log(model_file)
    testthat::expect_gt(nchar(parsed_model), 0L)

    if (!direct_ok) {
      outdir <- withr::local_tempdir()
      epa <- tc_run_epa_capture(
        epa_exec = epa_exec,
        ref_file = ref_file,
        query_file = query_file,
        tree_file = tree_file,
        model_arg = parsed_model,
        outdir = outdir
      )
      epa_log <- paste(c(epa$stdout, epa$stderr), collapse = "\n")
      testthat::expect_equal(
        epa$status,
        0L,
        info = paste("EPA parser fallback failed for", model_name)
      )
    }

    jplace <- file.path(outdir, "epa_result.jplace")
    testthat::expect_true(file.exists(jplace), info = model_name)
    testthat::expect_match(
      epa_log,
      "Specified model",
      info = paste("No model line for", model_name)
    )
    if (expected$has_base_freq) {
      testthat::expect_match(
        epa_log,
        "Base frequencies \\(empirical\\):",
        info = paste("Expected fitted base frequencies for", model_name)
      )
    }
    if (expected$has_rates) {
      testthat::expect_match(
        epa_log,
        "Substitution rates \\(ML\\):",
        info = paste("Expected fitted substitution rates for", model_name)
      )
    }
    if (expected$has_gamma) {
      testthat::expect_match(
        epa_log,
        "Rate heterogeneity:",
        info = paste("Expected fitted rate heterogeneity for", model_name)
      )
    }
  }
})

test_that("bayesant trains and predicts on generated fixtures", {
  tc_skip_if_no_pkg("BayesANT")
  tc_skip_if_no_pkg("Biostrings")
  x <- tc_generate_dataset(seed = 123L)
  tc_assert_dataset_shape(x)

  ref_file <- withr::local_tempfile(fileext = ".fasta")
  query_file <- withr::local_tempfile(fileext = ".fasta")
  tc_write_bayesant_reference(x$reference, ref_file)
  tc_write_plain_query(x$query, query_file)

  train <- BayesANT::read.BayesANT.data(
    ref_file,
    rank_names = optimotu.pipeline::tax_ranks()
  )
  model <- BayesANT::BayesANT(train, typeseq = "aligned", verbose = FALSE)

  out <- optimotu.pipeline::bayesant(
    query = query_file,
    model = model,
    ncpu = 1L,
    id_is_int = TRUE,
    n_top_taxa = 3L,
    min_prob = 0
  )

  checkmate::expect_names(
    names(out),
    must.include = c("seq_idx", "rank", "parent_taxonomy", "taxon", "prob")
  )
  testthat::expect_gt(nrow(out), 0L)
  testthat::expect_true(is.integer(out$seq_idx))
  testthat::expect_true(all(out$prob >= 0 & out$prob <= 1))
})
