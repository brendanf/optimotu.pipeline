# run_protax() is exercised with a no-op stub so we do not require a real
# ProtaxFungi install; the stub copies outdir/all.fa into the model directory
# for assertions. Use .cmd on Windows because system2()/CreateProcess cannot
# run shebang scripts there.
local_protax_stub <- function(envir = parent.frame()) {
  root <- withr::local_tempdir(
    pattern = "run_protax_test",
    .local_envir = envir
  )
  modeldir <- file.path(root, "model")
  dir.create(modeldir, recursive = TRUE)
  if (.Platform$OS.type == "windows") {
    script <- file.path(root, "runprotax_stub.cmd")
    writeLines(
      c(
        "@echo off",
        "copy /Y \"%~1\\all.fa\" \"%~2\\_captured_all.fa\" >NUL"
      ),
      script
    )
  } else {
    script <- file.path(root, "runprotax_stub.sh")
    writeLines(
      c(
        "#!/bin/sh",
        "set -e",
        "cp \"$1/all.fa\" \"$2/_captured_all.fa\""
      ),
      script
    )
    Sys.chmod(script, "0700")
  }
  list(
    root = root,
    modeldir = modeldir,
    script = script,
    captured = file.path(modeldir, "_captured_all.fa")
  )
}

test_that("run_protax returns character(0) for empty sequences", {
  stub <- local_protax_stub()
  outdir <- file.path(stub$root, "out_empty")
  bad_script <- file.path(stub$root, "never_run.sh")
  writeLines(c("#!/bin/sh", "exit 1"), bad_script)
  Sys.chmod(bad_script, "0700")
  got <- run_protax(
    Biostrings::DNAStringSet(),
    outdir = outdir,
    modeldir = stub$modeldir,
    ncpu = 1L,
    script = bad_script
  )
  expect_identical(got, character())
})
test_that("run_protax materializes DNAStringSet to all.fa for the stub", {
  stub <- local_protax_stub()
  outdir <- file.path(stub$root, "out_dss")
  seqs <- Biostrings::DNAStringSet(c(s1 = "ACGT", s2 = "NNNN"))
  run_protax(
    seqs,
    outdir = outdir,
    modeldir = stub$modeldir,
    ncpu = 1L,
    script = stub$script
  )
  expect_true(file.exists(stub$captured))
  got <- Biostrings::readDNAStringSet(stub$captured)
  expect_equal(as.character(got), as.character(seqs))
  expect_equal(names(got), names(seqs))
})
test_that("run_protax materializes a FASTA path to all.fa (not passthrough)", {
  stub <- local_protax_stub()
  outdir <- file.path(stub$root, "out_path")
  dir.create(outdir, recursive = TRUE)
  fa <- tempfile(fileext = ".fasta", tmpdir = stub$root)
  write_sequence(
    Biostrings::DNAStringSet(c(x = "AAAA", y = "TTTT")),
    fa
  )
  run_protax(
    fa,
    outdir = outdir,
    modeldir = stub$modeldir,
    ncpu = 1L,
    script = stub$script
  )
  got <- Biostrings::readDNAStringSet(stub$captured)
  expect_equal(names(got), c("x", "y"))
  expect_equal(unname(as.character(got)), c("AAAA", "TTTT"))
})
test_that("run_protax respects seq_idx for fastqindexr_index input", {
  stub <- local_protax_stub()
  outdir <- file.path(stub$root, "out_idx")
  fa <- tempfile(fileext = ".fasta", tmpdir = stub$root)
  write_sequence(
    Biostrings::DNAStringSet(c(p = "ACGT", q = "TGCA", r = "GGGG")),
    fa
  )
  idx <- fastqindexr::create_index(files = fa, type = "fasta")
  run_protax(
    idx,
    outdir = outdir,
    modeldir = stub$modeldir,
    ncpu = 1L,
    seq_idx = c(3L, 1L, 3L),
    script = stub$script
  )
  got <- Biostrings::readDNAStringSet(stub$captured)
  expect_equal(unname(as.character(got)), c("GGGG", "ACGT", "GGGG"))
  expect_equal(names(got), c("r", "p", "r"))
})
test_that("run_protax materializes named literal character sequences", {
  stub <- local_protax_stub()
  outdir <- file.path(stub$root, "out_char")
  run_protax(
    c(u = "AAA", v = "CCC"),
    outdir = outdir,
    modeldir = stub$modeldir,
    ncpu = 1L,
    script = stub$script
  )
  got <- Biostrings::readDNAStringSet(stub$captured)
  expect_equal(names(got), c("u", "v"))
  expect_equal(unname(as.character(got)), c("AAA", "CCC"))
})
