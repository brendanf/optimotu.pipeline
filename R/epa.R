find_epa_ng <- function() {
  find_executable("epa-ng")
}

#' Run EPA-ng
#' @param ref_msa ([`XStringSet`][Biostrings::XStringSet-class],
#' [`MultipleAlignment`][Biostrings::MultipleAlignment-class],
#' or `character`) reference alignment, or a file name
#' @param tree ([`phylo`][ape::read.tree] or `character`) reference tree, or
#' file name of a tree in Newick format
#' @param query ([`XStringSet`][Biostrings::XStringSet-class],
#' [`MultipleAlignment`][Biostrings::MultipleAlignment-class],
#' or `character`) query alignment, or a file name
#' @param outdir (`character`) output directory
#' @param model (`character`) substitution model specification, or a file name
#' containing the model specification.  EPA-ng can parse output files from
#' various programs, including IQ-TREE, RAxML, and FastTree.
#' @param ncpu (`integer`) number of threads to use
#' @param exec (`character`) path to the `epa-ng` executable
#' @param redo (`logical`) whether to overwrite existing output files
#' @param strip_inserts (`logical`) whether to strip insertions (lower-case
#' characters) from the query sequences before running EPA-ng
#' @return (`character`) path to the output `epa_result.jplace` file
#' @export
epa_ng <- function(
    ref_msa,
    tree,
    query,
    outdir = tempfile(),
    model,
    ncpu,
    exec = find_epa_ng(),
    redo = TRUE,
    strip_inserts = FALSE
) {
  checkmate::assert_flag(strip_inserts)

  ref_msa_file <- tempfile("reference", fileext = ".fasta")
  if (methods::is(ref_msa, "XStringSet")) {
    assertthat::assert_that(length(unique(Biostrings::width(ref_msa))) == 1)
  } else if (methods::is(ref_msa, "MultipleAlignment")) {
    ref_msa <- methods::as(ref_msa, "XStringSet")
  } else if (is.character(ref_msa)) {
    if (length(ref_msa) == 1 && file.exists(ref_msa)) {
      ref_msa_file <- ref_msa
    } else {
      assertthat::assert_that(length(unique(nchar(ref_msa))) == 1)
      if (has_alphabet(ref_msa, Biostrings::DNA_ALPHABET)) {
        ref_msa <- Biostrings::DNAStringSet(ref_msa)
      } else if (has_alphabet(ref_msa, Biostrings::RNA_ALPHABET)) {
        ref_msa <- Biostrings::RNAStringSet(ref_msa)
      } else if (has_alphabet(ref_msa, Biostrings::AA_ALPHABET)) {
        ref_msa <- Biostrings::AAStringSet(ref_msa)
      } else {
        stop("Unknown alphabet in reference alignment.")
      }
    }
  } else {
    stop("'ref_msa' should be an XStringSet, MultipleAlignment, character vector",
         "of aligned sequences, or filename.")
  }
  if (methods::is(ref_msa, "XStringSet")) {
    Biostrings::writeXStringSet(ref_msa, ref_msa_file)
    on.exit(file.remove(ref_msa_file))
  }

  query_file <- tempfile("query", fileext = ".fasta")
  if (methods::is(query, "XStringSet")) {
    assertthat::assert_that(length(unique(Biostrings::width(query))) == 1)
  } else if (methods::is(query, "MultipleAlignment")) {
    query <- methods::as(query, "XStringSet")
  } else if (is.character(query)) {
    if (length(query) == 1 && file.exists(query)) {
      query_file <- query
    } else {
      assertthat::assert_that(length(unique(nchar(query))) == 1)
      if (has_alphabet(query, Biostrings::DNA_ALPHABET)) {
        query <- Biostrings::DNAStringSet(query)
      } else if (has_alphabet(query, Biostrings::RNA_ALPHABET)) {
        query <- Biostrings::RNAStringSet(query)
      } else if (has_alphabet(query, Biostrings::AA_ALPHABET)) {
        query <- Biostrings::AAStringSet(query)
      } else {
        stop("Unknown alphabet in query alignment.")
      }
    }
  } else {
    stop("'query' should be an XStringSet, MultipleAlignment, character vector",
         "of aligned sequences, or filename.")
  }
  if (methods::is(query, "XStringSet")) {
    Biostrings::writeXStringSet(query, query_file)
  }

  is_gz <- endsWith(query_file, ".gz")
  if (is_gz | strip_inserts) {
    pre_file <- query_file
    query_file <- withr::local_tempfile(fileext=".fasta")
    precommand <- paste(if (is_gz) "zcat" else "cat", pre_file)
    if (strip_inserts) precommand <- paste(precommand, "| tr -d 'acgt'")
    precommand <- paste(precommand, ">", query_file)
    cat("running precommand:", precommand, "\n")
    pre_status <- system(precommand)
    stopifnot(pre_status == 0)
  }

  if (methods::is(tree, "phylo")) {
    tree_file <- tempfile("tree", tempdir(), ".tree")
    ape::write.tree(tree, tree_file)
    on.exit(file.remove(tree_file))
  } else if (is.character(tree) && file.exists(tree)) {
    tree_file <- tree
  } else {
    stop("'tree' should be a phylo object or a file name.")
  }

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  # if (length(model) == 1 && file.exists(model)) {
  #   model_file <- model
  # } else {
  #   model_file <- tempfile("info")
  #   on.exit(unlink(model_file))
  #   writeLines(model, model_file)
  # }

  args <- c(
    "--ref-msa", ref_msa_file,
    "--query", query_file,
    "--tree", tree_file,
    "--outdir", outdir,
    "--model", model
  )

  if (isTRUE(redo)) {
    args <- c(args, "--redo")
  }

  if (!missing(ncpu) && !is.null(ncpu)) {
    assertthat::assert_that(assertthat::is.count(ncpu))
    args <- c(args, "--threads", ncpu)
  }

  processx::run(exec, args = args, echo_cmd = TRUE, stdout = "")
  outfile <- file.path(outdir, "epa_result.jplace")
  checkmate::assert_file_exists(outfile, "r")
  outfile
}

#' Test whether a jplace object is valid
is_jplace <- function(jplace) {
  is.list(jplace) &&
    setequal(
      names(jplace),
      c("tree", "placements", "metadata", "version", "fields")
    )
}

#' Parse IQ-TREE model specification
#'
#' @details The current implementation can only parse a GTR+FU model.
#' @param model (`character`) model specification or a file name
parse_iqtree_log <- function(model) {
  checkmate::assert_character(model)
  if (checkmate::test_file_exists(model, access = "r")) {
    model = readLines(model)
  }
  model <- model[cumsum(grepl("FINALIZING TREE SEARCH", model)) == 1]
  rates <- model[grepl("Rate parameters:", model)]
  outmodel <- ""
  if (length(rates) == 1) {
    rates <- regmatches(rates, regexec("A-C: *([0-9.]+) +A-G: *([0-9.]+) +A-T: *([0-9.]+) +C-G: *([0-9.]+) +C-T: *([0-9.]+) +G-T: *([0-9.]+)", rates))[[1]][-1]
    rates <- as.numeric(rates)
    checkmate::assert_numeric(rates, lower = 0, finite = TRUE, len = 6, any.missing = FALSE)
    outmodel <- paste0(outmodel, "GTR{", paste0(rates, collapse = "/"), "}")
  }
  base_freq <- model[grepl("Base frequencies:", model)]
  if (length(base_freq) == 1) {
    base_freq <- regmatches(
      base_freq,
      regexec(
        pattern = "A: *([0-9.]+) +C: *([0-9.]+) +G: *([0-9.]+) +T: *([0-9.]+)",
        text = base_freq
      )
    )[[1]][-1]
    base_freq <- as.numeric(base_freq)
    checkmate::assert_numeric(base_freq, lower = 0, upper = 1, len = 4, any.missing = FALSE)
    outmodel <- paste0(outmodel, "+FU{", paste0(base_freq, collapse = "/"), "}")
  }
  alpha <- model[grepl("Gamma shape parameter:", model)]
  outmodel
}

parse_iqtree_model <- function(file = NULL, text = NULL) {
  if (is.null(text)) {
    checkmate::assert_file_exists(file, access = "r")
    text <- readLines(file)
  } else if (!is.null(file)) {
    stop("Cannot provide both 'file' and 'text'.")
  }
  checkmate::assert_character(text, any.missing = FALSE)
  if (length(text) == 1) text <- unlist(strsplit(text, "\n"))
  checkmate::assert_character(text[1], pattern = "^IQ-TREE")

  input_data_regex <-
    "^Input data: [0-9]+ sequences with [0-9]+ (amino-acid|nucleotide) sites$"
  input_data_line <- grep(input_line_regex, text, value = TRUE)
  checkmate::assert_character(input_data_line, len = 1)
  alphabet <- sub(input_data_regex, "\\1", input_data_line)

  model_regex <- "^Model of substitution: (.+)$"
  model_line <- grep(model_regex, text, value = TRUE)
  checkmate::assert_character(model_line, len = 1)
  model <- sub(model_regex, "\\1", model_line)

  model <- strsplit(model, split = "+", fixed = TRUE)


  rates <- text[grepl("Rate parameters:", text)]
  outmodel <- ""
  if (length(rates) == 1) {
    rates <- regmatches(rates, regexec("A-C: *([0-9.]+) +A-G: *([0-9.]+) +A-T: *([0-9.]+) +C-G: *([0-9.]+) +C-T: *([0-9.]+) +G-T: *([0-9.]+)", rates))[[1]][-1]
    rates <- as.numeric(rates)
    checkmate::assert_numeric(rates, lower = 0, finite = TRUE, len = 6, any.missing = FALSE)
    outmodel <- paste0(outmodel, "GTR{", paste0(rates, collapse = "/"), "}")
  }
  base_freq <- text[grepl("Base frequencies:", text)]
  if (length(base_freq) == 1) {
    base_freq <- regmatches(
      base_freq,
      regexec(
        pattern = "A: *([0-9.]+) +C: *([0-9.]+) +G: *([0-9.]+) +T: *([0-9.]+)",
        text = base_freq
      )
    )[[1]][-1]
    base_freq <- as.numeric(base_freq)
    checkmate::assert_numeric(base_freq, lower = 0, upper = 1, len = 4, any.missing = FALSE)
    outmodel <- paste0(outmodel, "+FU{", paste0(base_freq, collapse = "/"), "}")
  }
  alpha <- text[grepl("Gamma shape parameter:", text)]
  if (length(model) == 0) {
    stop("No model found in the log file.")
  }
  model
}

#' Run Gappa
#' @param jplace (`character`) path to the jplace file or a list containing
#' jplace data
#' @param taxonomy (`character`) path to the taxonomy file or a data frame
#' with two columns containing the taxonomy data; the first column contains
#' sequence IDs, while the second column contains comma-separated taxonomic
#' classifications.
#' @param outgroup (`character`) path to the outgroup file; or a character
#' vector containing the IDs of the outgroup sequences; or a character vector
#' containing the names of one or more taxa present in the taxonomy file.
#' @param ranks (`character`) vector of taxonomic ranks to use for the
#' assignment
#' @param ncpu (`integer`) number of threads to use
#' @param allow_file_overwriting (`logical`) whether to overwrite existing
#' output files
#' @param verbose (`logical`) whether to print verbose output
#' @param id_is_int (`logical`) whether the sequence IDs are integers
#' @return (`data.frame`) with columns:
#' - `seq_id` (character) the sequence ID
#' - `rank` (ordered factor) the taxonomic rank
#' - `parent_taxon` (character) the parent taxon
#' - `taxon` (character) the taxon
#' - `prob` (numeric) the probability of the assignment
#' @export
gappa_assign <- function(
    jplace,
    taxonomy,
    outgroup,
    ranks = unknown_ranks(),
    ncpu = NULL,
    allow_file_overwriting = TRUE,
    verbose = FALSE,
    id_is_int = FALSE
) {
  gappa <- find_executable("gappa")
  args <- c("examine", "assign", "--per-query-results")

  if (checkmate::test_file_exists(jplace, access = "r")) {
    jplace_file <- jplace
  } else if (is_jplace(jplace)) {
    jplace_file <- withr::local_tempfile(fileext = ".jplace")
    jsonlite::write_json(jplace, jplace_file, auto_unbox = TRUE)
  } else {
    stop("'jplace' should be a filename or a list containing jplace data")
  }
  args <- c(args, "--jplace-path", jplace_file)

  if (checkmate::test_file_exists(taxonomy, access = "r")) {
    taxonomy_file <- taxonomy
  } else {
    checkmate::assert_data_frame(taxonomy, min.cols = 2, max.cols = 2)
    taxonomy_file <- withr::local_tempfile(fileext = ".txt")
    readr::write_tsv(taxonomy, taxonomy_file, col_names = FALSE)
  }
  args <- c(args, "--taxon-file", taxonomy_file)

  if (checkmate::test_file_exists(outgroup, access = "r")) {
    outgroup_file <- outgroup
  } else {
    checkmate::assert_character(outgroup, any.missing = FALSE, min.chars = 1,
                                 min.len = 1)
    if (!is.data.frame(taxonomy)) {
      taxonomy <- readr::read_tsv(taxonomy_file, col_names = FALSE, col_types = "cc")
    }
    if (!all(outgroup) %in% taxonomy[[1]]) {
      outgroup <- taxonomy[[1]][grepl(paste(outgroup, collapse = "|"), taxonomy[[2]])]
      if (length(outgroup) == 0) {
        stop("Outgroup not found in the taxonomy file.")
      }
    }
    outgroup_file <- withr::local_tempfile(fileext = ".txt")
    writeLines(outgroup, outgroup_file)
  }
  args <- c(args, "--root-outgroup", outgroup_file)

  checkmate::assert_character(ranks, any.missing = FALSE)
  args <- c(args, "--ranks-string", paste(ranks, collapse = "|"))

  checkmate::assert_flag(allow_file_overwriting, null.ok = TRUE)
  if (isTRUE(allow_file_overwriting)) args <- c(args, "--allow-file-overwriting")

  checkmate::assert_flag(verbose, null.ok = TRUE)
  if (isTRUE(verbose)) args <- c(args, "--verbose")

  checkmate::assert_flag(id_is_int)

  checkmate::assert_count(ncpu, null.ok = TRUE)
  if (!is.null(ncpu)) args <- c(args, "--threads", ncpu)

  out_file <- withr::local_tempfile(fileext = "_per_query.tsv")
  checkmate::assert_path_for_output(out_file)
  args <- c(
    args,
    "--out-dir", dirname(out_file),
    "--file-prefix", sub("per_query.tsv$", "", basename(out_file))
  )

  processx::run(gappa, args = args, echo_cmd = TRUE, stdout = "")
  parse_gappa_per_query(out_file, ranks, id_is_int)
}

#' Parse Gappa per-query results
#' @param per_query (`character`) path to the per-query results file
#' @param ranks (`character`) vector of taxonomic ranks
#' @param id_is_int (`logical`) whether the sequence IDs are integers
#' @keywords internal
parse_gappa_per_query <- function(per_query, ranks, id_is_int = FALSE) {
  checkmate::assert_file(per_query, "r")
  checkmate::assert_character(ranks, any.missing = FALSE)
  checkmate::assert_flag(id_is_int)
  pq <- readr::read_tsv(per_query, col_types = if (id_is_int) "innnnc" else "cnnnnc") |>
    dplyr::mutate(
      rank = factor(stringr::str_count(taxopath, ";") + 1, levels = seq_along(ranks), labels = ranks, ordered = TRUE),
      taxopath = chartr(";", ",", taxopath)
    )
  pq <- dplyr::bind_rows(
    dplyr::filter(pq, fract > 0, rank != dplyr::last(ranks)) |>
      dplyr::transmute(
        name,
        rank = factor(ranks[as.integer(rank) + 1L], levels = ranks, ordered = TRUE),
        parent_taxon = taxopath,
        taxon = NA_character_,
        prob = fract
      ) |>
      dplyr::filter(prob > 0),
    tidyr::extract(pq, taxopath, c("parent_taxon", "taxon"), , regex = "(?:(.+),)?([^,]+)$") |>
      dplyr::transmute(
        name,
        rank,
        parent_taxon = dplyr::na_if(parent_taxon, ""),
        taxon,
        prob = afract
      )
  ) |>
    dplyr::arrange(name, rank)

  names(pq)[1] <- if (id_is_int) "seq_idx" else "seq_id"
  pq
}
