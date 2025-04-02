#' Find OTUs which contain sequences with _any_ probability of being a target taxon
#' @param target_taxa (`character`) target taxa
#' @param asv_all_tax_prob (`data.frame`) ASV taxonomic probabilities, as
#' returned by `parse_protax_nameprob()` or `run_protax_animal()`. Should
#' include columns `seq_id`, `rank`, `taxon`, and `prob`.
#' @param asv_taxonomy (`data.frame`) ASV taxonomy; as `asv_all_tax_prob`, but
#' with only one taxon per rank per ASV, and `prob` is not required.
#' @param otu_taxonomy (`data.frame`) OTU taxonomy; as `asv_taxonomy`, but for OTUs.
#' @return `data.frame` with columns `seq_id`, `asv_seq_id`, `otu_taxon`, `rank`,
#' `protax_taxon`, and `protax_prob`.
#' @export
find_target_taxa <- function(target_taxa, asv_all_tax_prob, asv_taxonomy, otu_taxonomy) {
  # avoid R CMD check NOTE: no visible binding for global variable
  seq_id <- otu_taxon <- taxon <- prob <- protax_taxon <- asv_seq_id <- NULL

  checkmate::assert_character(target_taxa)
  checkmate::assert_data_frame(asv_all_tax_prob)
  checkmate::assert_subset(c("seq_id", "rank", "taxon", "prob"), colnames(asv_all_tax_prob))
  checkmate::assert_data_frame(asv_taxonomy)
  checkmate::assert_subset(c("seq_id", "rank", "taxon"), colnames(asv_taxonomy))
  checkmate::assert_data_frame(otu_taxonomy)
  checkmate::assert_subset(c("seq_id", "rank", "taxon"), colnames(otu_taxonomy))

  asv_otu_key <-
    dplyr::inner_join(
      dplyr::select(asv_taxonomy, asv_seq_id = seq_id, !!tip_rank_var()),
      dplyr::select(otu_taxonomy, seq_id, !!tip_rank_var()),
      by = tip_rank()
    ) |>
    dplyr::select(-(!!tip_rank_var()))

  otu_long_taxonomy <- tidyr::pivot_longer(
    otu_taxonomy,
    all_of(tax_ranks()),
    names_to = "rank",
    values_to = "otu_taxon",
    names_transform = list(rank = rank2factor)
  ) |>
    dplyr::select(seq_id, rank, otu_taxon)

  dplyr::select(
    asv_all_tax_prob,
    asv_seq_id = seq_id,
    rank,
    protax_taxon = taxon,
    protax_prob = prob
  ) |>
    dplyr::inner_join(asv_otu_key, by = "asv_seq_id") |>
    dplyr::group_by(seq_id) |>
    dplyr::filter(any(protax_taxon %in% target_taxa)) |>
    dplyr::left_join(otu_long_taxonomy, by = c("seq_id", "rank")) |>
    dplyr::ungroup() |>
    dplyr::arrange(seq_id, asv_seq_id, rank) |>
    dplyr::select(seq_id, asv_seq_id, otu_taxon, rank, everything())
}
