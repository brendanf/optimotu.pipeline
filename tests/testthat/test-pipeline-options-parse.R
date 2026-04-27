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

  fasta <- withr::local_tempfile(fileext = ".fa")
  table <- withr::local_tempfile(fileext = ".xlsx")
  writeLines(">x\nACGT", fasta)
  writeLines("placeholder", table)

  expect_error(
    optimotu.pipeline:::parse_added_reference(
      list(added_reference = list(fasta = fasta))
    ),
    "both must be given"
  )

  optimotu.pipeline:::parse_added_reference(
    list(added_reference = list(list(fasta = fasta), list(table = table)))
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
