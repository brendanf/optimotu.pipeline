test_that("consensus_columns works on legacy list MultipleAlignment input", {
  aln <- Biostrings::DNAMultipleAlignment(
    c(
      seq1 = "ACGT....ACGT",
      seq2 = "ACGA....TCGT"
    )
  )
  rf <- Biostrings::BString("xxxx....xxxx")
  out <- consensus_columns(list(alignment = aln, GC = list(RF = rf)))
  expect_s4_class(out, "DNAStringSet")
  expect_equal(unname(Biostrings::width(out)), c(8L, 8L))
  expect_equal(as.character(out[["seq1"]]), "ACGTACGT")
  expect_equal(as.character(out[["seq2"]]), "ACGATCGT")
})

test_that("consensus_columns handles RF with no insert columns", {
  aln <- Biostrings::DNAMultipleAlignment(
    c(
      seq1 = "ACGTACGT",
      seq2 = "ACGATCGT"
    )
  )
  rf <- Biostrings::BString("xxxxxxxx")
  out <- consensus_columns(list(alignment = aln, GC = list(RF = rf)))
  expect_s4_class(out, "DNAStringSet")
  expect_equal(as.character(out[["seq1"]]), "ACGTACGT")
  expect_equal(as.character(out[["seq2"]]), "ACGATCGT")
})

test_that("consensus_columns works on StockholmDNAMultipleAlignment when available", {
  skip_if_not_installed("inferrnal")
  aln <- inferrnal::StockholmDNAMultipleAlignment(
    c(
      seq1 = "ACGT....ACGT",
      seq2 = "ACGA....TCGT"
    ),
    GC = c(RF = "xxxx....xxxx")
  )
  out <- consensus_columns(aln)
  expect_s4_class(out, "DNAStringSet")
  expect_equal(unname(Biostrings::width(out)), c(8L, 8L))
  expect_equal(as.character(out[["seq1"]]), "ACGTACGT")
  expect_equal(as.character(out[["seq2"]]), "ACGATCGT")
})

test_that("consensus_columns converts StockholmRNA to DNAStringSet", {
  skip_if_not_installed("inferrnal")
  aln <- inferrnal::StockholmRNAMultipleAlignment(
    c(
      seq1 = "ACGU....ACGU",
      seq2 = "ACGA....UCGU"
    ),
    GC = c(RF = "xxxx....xxxx")
  )
  out <- consensus_columns(aln)
  expect_s4_class(out, "DNAStringSet")
  expect_equal(as.character(out[["seq1"]]), "ACGTACGT")
  expect_equal(as.character(out[["seq2"]]), "ACGATCGT")
})
