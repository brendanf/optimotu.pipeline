# Deterministic synthetic fixtures for classifier integration tests.

tc_bases <- c("A", "C", "G", "T")

tc_random_dna <- function(n) {
  paste(sample(tc_bases, n, replace = TRUE), collapse = "")
}

tc_mutate_sequence <- function(seq, frac) {
  chars <- strsplit(seq, "", fixed = TRUE)[[1]]
  n_mut <- max(1L, as.integer(round(length(chars) * frac)))
  idx <- sample.int(length(chars), n_mut, replace = FALSE)
  for (i in idx) {
    chars[[i]] <- sample(tc_bases[tc_bases != chars[[i]]], 1L)
  }
  paste(chars, collapse = "")
}

tc_generate_dataset <- function(seed = 20260423L) {
  withr::with_seed(seed, {
    ranks <- optimotu.pipeline::tax_ranks()
    id <- 1L
    root <- tc_random_dna(300L)
    rows <- vector("list", 81L)
    row_i <- 1L

    for (f in seq_len(3L)) {
      family_seq <- tc_mutate_sequence(root, 0.15)
      family <- sprintf("Family_%02d", f)
      for (g in seq_len(3L)) {
        genus_seq <- tc_mutate_sequence(family_seq, 0.10)
        genus <- sprintf("Genus_%02d_%02d", f, g)
        for (s in seq_len(3L)) {
          species_seq <- tc_mutate_sequence(genus_seq, 0.05)
          species <- sprintf("Species_%02d_%02d_%02d", f, g, s)
          for (i in seq_len(3L)) {
            individual <- tc_mutate_sequence(species_seq, 0.01)
            rows[[row_i]] <- tibble::tibble(
              seq_id = as.character(id),
              seq = individual,
              kingdom = "Life",
              phylum = "P_Test",
              class = "C_Test",
              order = "O_Test",
              family = family,
              genus = genus,
              species = species,
              set = if (i <= 2L) "reference" else "query"
            )
            row_i <- row_i + 1L
            id <- id + 1L
          }
        }
      }
    }

    all <- dplyr::bind_rows(rows)
    all <- all[, c("seq_id", "seq", ranks, "set"), drop = FALSE]
    list(
      root = root,
      all = all,
      reference = all[all$set == "reference", , drop = FALSE],
      query = all[all$set == "query", , drop = FALSE]
    )
  })
}

tc_assert_dataset_shape <- function(x) {
  testthat::expect_equal(nrow(x$reference), 54L)
  testthat::expect_equal(nrow(x$query), 27L)
  testthat::expect_equal(dplyr::n_distinct(x$all$family), 3L)
  testthat::expect_equal(dplyr::n_distinct(x$all$genus), 9L)
  testthat::expect_equal(dplyr::n_distinct(x$all$species), 27L)
}

tc_write_fasta <- function(seq, id, file) {
  dss <- Biostrings::DNAStringSet(stats::setNames(seq, id))
  Biostrings::writeXStringSet(dss, file)
  file
}

tc_sintax_header <- function(df) {
  paste0(
    df$seq_id,
    ";tax=k:", df$kingdom,
    ",p:", df$phylum,
    ",c:", df$class,
    ",o:", df$order,
    ",f:", df$family,
    ",g:", df$genus,
    ",s:", df$species
  )
}

tc_bayesant_header <- function(df) {
  paste(
    df$seq_id,
    paste(
      c(
        "Root",
        df$kingdom,
        df$phylum,
        df$class,
        df$order,
        df$family,
        df$genus,
        df$species
      ),
      collapse = ";"
    )
  )
}

tc_taxonomy_path <- function(df) {
  apply(
    df[, c("kingdom", "phylum", "class", "order", "family", "genus", "species")],
    1L,
    paste,
    collapse = ","
  )
}

tc_write_sintax_reference <- function(df, file) {
  tc_write_fasta(df$seq, tc_sintax_header(df), file)
}

tc_write_plain_query <- function(df, file) {
  tc_write_fasta(df$seq, df$seq_id, file)
}

tc_write_bayesant_reference <- function(df, file) {
  tc_write_fasta(df$seq, tc_bayesant_header(df), file)
}

tc_write_gappa_taxonomy <- function(df, file) {
  readr::write_tsv(
    tibble::tibble(seq_id = df$seq_id, taxonomy = tc_taxonomy_path(df)),
    file,
    col_names = FALSE
  )
  file
}

tc_iqtree_model_matrix <- function() {
  c("JC", "HKY", "GTR", "GTR+F", "HKY+F", "GTR+G", "GTR+I+G")
}

tc_parse_iqtree_expected <- function(iqtree_file) {
  txt <- readLines(iqtree_file, warn = FALSE)
  model_line <- grep("^Model of substitution: ", txt, value = TRUE)
  base_freq <- grep("^Base frequencies:", txt, value = TRUE)
  rates <- grep("^Rate parameters:", txt, value = TRUE)
  alpha <- grep("^Gamma shape parameter:", txt, value = TRUE)
  list(
    model_line = if (length(model_line) > 0) model_line[[1]] else "",
    has_base_freq = length(base_freq) > 0,
    has_rates = length(rates) > 0,
    has_gamma = length(alpha) > 0
  )
}

tc_run_epa_capture <- function(
    epa_exec,
    ref_file,
    query_file,
    tree_file,
    model_arg,
    outdir
) {
  processx::run(
    command = epa_exec,
    args = c(
      "--ref-msa",
      ref_file,
      "--query",
      query_file,
      "--tree",
      tree_file,
      "--outdir",
      outdir,
      "--model",
      model_arg,
      "--redo",
      "--threads",
      "1"
    ),
    stdout = "|",
    stderr = "|",
    error_on_status = FALSE
  )
}

tc_find_exec <- function(candidates) {
  for (cmd in candidates) {
    path <- Sys.which(cmd)
    if (nzchar(path)) {
      return(path)
    }
  }
  ""
}

tc_find_gappa <- function() {
  sys <- "/usr/local/bin/gappa"
  if (file.exists(sys) && file.access(sys, 1) == 0) {
    return(sys)
  }
  tc_find_exec("gappa")
}

tc_skip_if_no_exec <- function(candidates, name = NULL) {
  path <- tc_find_exec(candidates)
  if (!nzchar(path)) {
    label <- name %||% paste(candidates, collapse = " or ")
    testthat::skip(paste0("Missing executable: ", label))
  }
  path
}

tc_skip_if_no_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    testthat::skip(paste0("Missing package: ", pkg))
  }
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
