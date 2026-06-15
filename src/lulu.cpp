// implementation of LULU algorithm for a long OTU table
// Written without reference to the LULU code, based only on algorithm
// description in the manuscript and documentation.

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include "accessor.h"


// struct to store match information about a pair of sequences
// abund_ratio is either the sum of abundance ratios in all samples where the
//  sequences co-occur (if use_mean_abundance = true) or the minimum abundance
//  ratio (if use_mean_abundance_ratio = false).
// nboth is the number of samples in which they co-occur.
// It is assumed that the first sequence is the potential child, i.e. the one
// which is less prevalent, or if tied then less abundant.
struct match_info {
  double abund_ratio = 0;
  int nboth = 0;

  match_info(int nread1, int nread2) :
    abund_ratio(double(nread2) / double(nread1)),
    nboth(1) {}

  match_info() {}

  void add_match(int nread1, int nread2, bool use_mean_abundance_ratio) {
    double new_abund_ratio = double(nread2) / double(nread1);
    if (nboth == 0) {
      abund_ratio = new_abund_ratio;
    } else if (use_mean_abundance_ratio) {
      abund_ratio += new_abund_ratio;
    } else {
      abund_ratio = std::min(abund_ratio, new_abund_ratio);
    }
    ++nboth;
  }
};

struct match_info_data {
  std::map<std::pair<int, int>, match_info> data;
  const bool use_mean;

  match_info_data(bool use_mean_abundance_ratio) :
    use_mean(use_mean_abundance_ratio) {}


  void add_match(int id1, int id2, int nread1, int nread2) {
    // ensure that id1 is the larger ID (i.e., the potential child)
    // this means that we only need to check one direction for matches
    if (id2 > id1) {
      std::swap(id1, id2);
      std::swap(nread1, nread2);
    }
    data[std::make_pair(id1, id2)].add_match(nread1, nread2, use_mean);
  }
};

// overload for core implementation
Rcpp::DataFrame lulu_map_impl(
    Rcpp::IntegerVector seq_idx_out,
    match_info_data & match_info,
    std::vector<int> & total_occurrence,
    double min_abundance_ratio = 1.0,
    double min_cooccurrence_ratio = 0.95,
    bool use_mean_abundance_ratio = false,
    int verbose = 0
) {

  std::vector<int> lulu_map(total_occurrence.size(), NA_INTEGER);
  for (int i : seq_idx_out) {
    lulu_map[i] = i;
  }

  for (const auto & mi : match_info.data) {
    if (verbose > 0) {
      Rcpp::Rcerr << "Considering match pair (" << mi.first.first
                  << ", " << mi.first.second << ") with "
                  << total_occurrence[mi.first.first] << " and "
                  << total_occurrence[mi.first.second] << " occurrences and "
                  << mi.second.nboth << " co-occurrences"
                  << std::endl;
    }
    // if the potential child has already been denoised, skip
    if (lulu_map[mi.first.first] != mi.first.first) {
      if (verbose > 0) {
        Rcpp::Rcerr << "seq " << mi.first.first
                    << " already mapped to seq " << lulu_map[mi.first.first]
                    << "; skipping" << std::endl;
      }
      continue;
    }

    // check the co-occurrence ratio
    if (mi.second.nboth < min_cooccurrence_ratio * total_occurrence[mi.first.first]) {
      if (verbose > 0) {
        Rcpp::Rcerr << "co-occurrence ratio " << mi.second.nboth
                    << " / " << total_occurrence[mi.first.first]
                    << " = " << double(mi.second.nboth) / double(total_occurrence[mi.first.first])
                    << " is less than minimum " << min_cooccurrence_ratio
                    << "; skipping" << std::endl;
      }
      continue;
    } else if (verbose > 1) {
      Rcpp::Rcerr << "co-occurrence ratio " << mi.second.nboth
                  << " / " << total_occurrence[mi.first.first]
                  << " = " << double(mi.second.nboth) / double(total_occurrence[mi.first.first])
                  << " is greater than or equal to minimum " << min_cooccurrence_ratio
                  << std::endl;
    }

    // check the abundance ratio
    double abundance_ratio = mi.second.abund_ratio;
    if (use_mean_abundance_ratio) {
      // the object has accumulated the sum, so we need to divide.
      abundance_ratio /= mi.second.nboth;
    }
    if (abundance_ratio > min_abundance_ratio) {
      if (verbose > 1) {
        Rcpp::Rcerr << (use_mean_abundance_ratio ? "mean abundance ratio " : "abundance ratio ")
                    << abundance_ratio << " greater than minimum "
                    << min_abundance_ratio << std::endl;
      }
      if (verbose > 0) {
        Rcpp::Rcerr << "Mapping child " << mi.first.first
                    << " to parent " << mi.first.second
                    << std::endl;
      }
      lulu_map[mi.first.first] = mi.first.second;
    } else {
      if (verbose > 0) {
        Rcpp::Rcerr << (use_mean_abundance_ratio ? "mean abundance ratio " : "abundance ratio ")
                    << abundance_ratio << " less than or equal to minimum " << min_abundance_ratio
                    << "; skipping" << std::endl;
      }
    }
  }

  Rcpp::IntegerVector lulu_idx_out(seq_idx_out.size());

  for (R_xlen_t i = 0; i < seq_idx_out.size(); i++) {
    int j = seq_idx_out[i];
    while (lulu_map[j] != j) {
      j = lulu_map[j];
    }
    lulu_idx_out[i] = j;
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("seq_idx") = seq_idx_out,
    Rcpp::Named("lulu_idx") = lulu_idx_out
  );
}

//' LULU secondary denoising
//'
//' The "match" inputs are intended to be calculated by pairwise distance
//' calculations within each sample, the results of which are concatenated
//' together. This requires some duplicate calculations for pairs that occur in
//' many samples, but prevents the need to calculate distances for pairs that
//' never co-occur, which is a large number of pairs.
//'
//' @param match_id1 (`integer`) vector of sequence IDs for the first sequence
//' in each match; sequence IDs should be assigned such that when the sequences
//' are in ID order, they are sorted by decreasing occurrence, with ties broken
//' by decreasing abundance.
//' @param match_id2 (`integer`) vector of sequence IDs for the second sequence
//' in each match
//' @param match_nread1 (`integer`) vector of read counts for the first sequence
//' in each match
//' @param match_nread2 (`integer`) vector of read counts for the second
//' sequence in each match
//' @param match_dist (`numeric`) vector of distances between the two sequences
//' in each match
//' @param seq_idx (`integer`) vector of sequence IDs. Each ID should occur
//' once per sample where the sequence occurs.
//' @param nread (`integer`) vector of read counts for each sequence in each
//' sample
//' @param max_dist (`numeric`) maximum distance for two sequences to be
//' considered as parent-child
//' @param min_abundance_ratio (`numeric`) minimum abundance ratio for two
//' sequences to be considered as parent-child
//' @param min_cooccurrence_ratio (`numeric`) minimum co-occurrence ratio for
//' two sequences to be considered
//' @param use_mean_abundance_ratio (`logical`) if `TRUE`, `min_abundance_ratio`
//' is interpreted as a minimum value for the mean of the abundance ratios in
//' all samples where the two sequences co-occur.  Otherwise (default) it is
//' interpreted as a minimum value for the the abundance ratio of all samples
//' where the two sequences co-occur.
//' @param verbose (`integer`) level of verbosity. At level 0 (default),
//' no messages are printed.
//' @returns a two-column `data.frame` with columns `seq_idx` and `lulu_idx`.
//' `seq_idx` includes all values which occur in the `seq_idx` argument, and
//' `lulu_idx` gives the index of the denoised sequence.
//'
// [[Rcpp::export]]
Rcpp::DataFrame lulu_map_impl(
  Rcpp::IntegerVector match_id1,
  Rcpp::IntegerVector match_id2,
  Rcpp::IntegerVector match_nread1,
  Rcpp::IntegerVector match_nread2,
  Rcpp::NumericVector match_dist,
  Rcpp::IntegerVector seq_idx,
  Rcpp::IntegerVector nread,
  double max_dist,
  double min_abundance_ratio = 1.0,
  double min_cooccurrence_ratio = 0.95,
  bool use_mean_abundance_ratio = false,
  int verbose = 0
) {
  // check that all match_* vectors are the same length
  if (match_id1.size() != match_id2.size() ||
      match_id1.size() != match_nread1.size() ||
      match_id1.size() != match_nread2.size() ||
      match_id1.size() != match_dist.size()) {
    Rcpp::stop("All match_* vectors must be the same length");
  }
  //check that seq_idx and seq_nsample are the same length
  if (seq_idx.size() != nread.size()) {
    Rcpp::stop("seq_idx and nread must be the same length");
  }

  // check that the sequences are correctly sorted by occurrence and abundance
  std::vector<int> total_occurrence;
  std::vector<int> total_abundance;
  int n_seq_idx = 0;
  for (int i = 0; i < seq_idx.size(); i++) {
    if (seq_idx[i] >= (int)total_occurrence.size()) {
      total_occurrence.resize(seq_idx[i] + 1, 0);
      total_abundance.resize(seq_idx[i] + 1, 0);
    }
    if (total_occurrence[seq_idx[i]] == 0) {
      n_seq_idx++;
    }
    total_occurrence[seq_idx[i]]++;
    total_abundance[seq_idx[i]] += nread[i];
  }

  Rcpp::IntegerVector seq_idx_out(n_seq_idx);
  Rcpp::IntegerVector lulu_idx_out(n_seq_idx);
  int j = 0;

  for (std::size_t i = 0; i < total_occurrence.size(); i++) {
    if (total_occurrence[i] > 0) {
      seq_idx_out[j] = i;
      j++;
    }
  }

  for (R_xlen_t i = 1; i < seq_idx_out.size(); i++) {
    if (total_occurrence[seq_idx_out[i]] > total_occurrence[seq_idx_out[i - 1]]) {
      Rcpp::stop(
        "seq_idx %d has %d occurences, greater than seq_idx %d with %d.",
        seq_idx_out[i],
        total_occurrence[seq_idx_out[i]],
        seq_idx_out[i - 1],
        total_occurrence[seq_idx_out[i - 1]]
      );
    }
    if (total_occurrence[seq_idx_out[i]] == total_occurrence[seq_idx_out[i - 1]] &&
        total_abundance[seq_idx_out[i]] > total_abundance[seq_idx_out[i - 1]]) {
      Rcpp::stop("seq_idx must be sorted by decreasing occurrence, with ties"
                 "broken by decreasing abundance");
    }
  }


  match_info_data match_info(use_mean_abundance_ratio);
  for (int i = 0; i < match_id1.size(); i++) {
    if (match_dist[i] <= max_dist) {
      match_info.add_match(match_id1[i], match_id2[i], match_nread1[i], match_nread2[i]);
    }
  }

  return lulu_map_impl(
    seq_idx_out,
    match_info,
    total_occurrence,
    min_abundance_ratio,
    min_cooccurrence_ratio,
    use_mean_abundance_ratio,
    verbose
  );
}

Rcpp::RObject tar_read(Rcpp::String name) {
  Rcpp::Environment targets = Rcpp::Environment::namespace_env("targets");
  Rcpp::Function tar_runtime_object = targets["tar_runtime_object"];
  Rcpp::Environment tar_runtime = tar_runtime_object();
  Rcpp::RObject meta_raw = tar_runtime["meta"];
  if (meta_raw.isNULL()) {
    // targets pipeline is not running, so we are allowed to use tar_read
    Rcpp::Function tar_read_raw = targets["tar_read_raw"];
    return tar_read_raw(name);
  } else {
    Rcpp::Environment meta = Rcpp::as<Rcpp::Environment>(meta_raw);
    Rcpp::Function exists_record = meta["exists_record"];
    Rcpp::LogicalVector exists = exists_record(name);
    if (exists[0] == FALSE) return R_NilValue;
    Rcpp::Function get_record = meta["get_record"];
    Rcpp::Environment record = get_record(name);
    Rcpp::Function record_bootstrap_store = targets["record_bootstrap_store"];
    Rcpp::RObject store = record_bootstrap_store(record);
    Rcpp::Function record_bootstrap_file = targets["record_bootstrap_file"];
    Rcpp::RObject file = record_bootstrap_file(record);
    Rcpp::Function store_read_object = targets["store_read_object"];
    return store_read_object(store, file);
  }
}


//' LULU secondary denoising for "big" data targets pipeline
//'
//' This version of `lulu_map()` is intended for integration in a `targets`
//' pipeline, where it should be used inside a target with option
//' `retrieval = "none"`. The assumption is that the LULU match table and/or the
//' OTU occurrence table are split into multiple files on disk using
//' dynamic and/or static branching in targets, and that it may be too memory
//' intensive to simply load all of the files to generate the full tables in
//' memory. Instead, they are processed one at a time, with calls to `gc()`
//' after each one. This saves memory in two ways: first, it avoids holding the
//' full data in both its R representation and in its internal C++
//' simultaneously; second, the internal C++ structures are smaller
//' than the R tables.
//'
//' @param match_table_names (`character` vector) fully resolved target names
//'   for the lulu match table, as generated by `lulu_distmx()`.
//' @param otu_table_names (`character` vector) fully resolved target names for
//'   the OTU occurrence table.
//' @inheritParams lulu_map_impl
// [[Rcpp::export]]
Rcpp::DataFrame lulu_map_lowmem_impl(
  Rcpp::CharacterVector otu_table_names,
  Rcpp::CharacterVector match_table_names,
  double max_dist,
  double min_abundance_ratio = 1.0,
  double min_cooccurrence_ratio = 0.95,
  bool use_mean_abundance_ratio = false,
  int verbose = 0

) {
  // First read the OTU table(s) to get totals
  std::vector<int> total_occurrences;
  std::vector<size_t> total_abundance;

  for (R_xlen_t i = 0; i < otu_table_names.size(); ++i) {
    Rcpp::String otu_table_name(otu_table_names[i]);
    if (verbose) {
      Rcpp::Rcerr << "Reading OTU table " << otu_table_name.get_cstring()
                  << "\n  Collecting garbage..." << std::flush;
    }
    R_gc();
    if (verbose) {
      Rcpp::Rcerr << "done.\n  Counting occurrences..." << std::flush;
    }
    Rcpp::RObject otu_table = tar_read(otu_table_name);
    Rcpp::IntegerVector seq_idx =
      integer_column(otu_table, "seq_idx", otu_table_name.get_cstring());
    Rcpp::IntegerVector nread =
      integer_column(otu_table, "nread", otu_table_name.get_cstring());

    for (R_xlen_t j = 0; j < seq_idx.size(); ++j) {
      if (seq_idx[j] >= (R_xlen_t)total_occurrences.size()) {
        total_occurrences.resize(seq_idx[j] + 1, 0);
        total_abundance.resize(seq_idx[j] + 1, 0);
      }
      total_occurrences[seq_idx[j]]++;
      total_abundance[seq_idx[j]] += nread[j];
    }
    if (verbose) Rcpp::Rcerr << "done." << std::endl;
  }
  if (verbose) Rcpp::Rcerr << "Collecting garbage..." << std::flush;
  R_gc();
  if (verbose) {
    Rcpp::Rcerr << "done.\nCounting used seq_idx values..." << std::flush;
  }
  // Count the seq indices that actually occur.
  std::size_t n_seq_idx = 0;
  for (int n : total_occurrences) {
    if (n > 0) ++n_seq_idx;
  }
  if (verbose) {
    Rcpp::Rcerr << "done.\nInitializing reverse map..." << std::flush;
  }
  // Create forward and reverse maps for sequence indices.
  // fwd_map[seq_idx] gives the position that a given seq_idx has in the
  //   ordering
  // rev_map[i] gives the seq_idx correcponding to a position in the ordering
  // Both are 0-indexed, but in typical usage (seq_idx are 1-based from R),
  // the maps will account for this because seq_idx = 0 is skipped because it
  // doesn't occur.

  // fwd_map is only used in C++
  std::vector<int> fwd_map(total_occurrences.size(), NA_INTEGER);

  // rev_map, nonempty_*, and order are used in calls to R API functions
  Rcpp::IntegerVector rev_map(n_seq_idx);
  Rcpp::IntegerVector nonempty_occurrences(n_seq_idx);
  Rcpp::IntegerVector nonempty_abundance(n_seq_idx);
  Rcpp::IntegerVector order(n_seq_idx);

  // i indexes over the
  int i = 0, j = 0;
  for (int n : total_occurrences) {
    if (n > 0) {
      rev_map[i] = j;
      nonempty_occurrences[i] = n;
      nonempty_abundance[i] = total_abundance[j];
      ++i;
    }
    ++j;
  }
  if (verbose) {
    Rcpp::Rcerr << "done.\nSorting reverse map..." << std::flush;
  }

  // C API for R's base::order()
  // It uses radix sort for integers so it is MUCH faster than a naive C++
  // implementation using std::sort with a custom comparator.
  R_orderVector(
    INTEGER(order),
    n_seq_idx,
    Rf_lang2(nonempty_occurrences, nonempty_abundance),
    FALSE,
    TRUE
  );

  rev_map = rev_map[order];

  if (verbose) {
    Rcpp::Rcerr << "done.\nInverting to form forward map..." << std::flush;
  }

  // invert the rev_map to get the fwd_map
  i = 0;
  for (int seq_idx : rev_map) {
    fwd_map[seq_idx] = i++;
  }

  if (verbose) Rcpp::Rcerr << "done." << std::endl;

  // Now read the match table(s) to count co-occurrences and relative abundances
  match_info_data mid(use_mean_abundance_ratio);

  for (R_xlen_t i = 0; i < match_table_names.size(); ++i) {
    Rcpp::String match_table_name(match_table_names[i]);
    if (verbose) {
      Rcpp::Rcerr << "Reading match table " << match_table_name.get_cstring()
                  << "\n  Collecting garbage.." << std::flush;
    }
    R_gc();
    if (verbose) {
      Rcpp::Rcerr << "done.\n  Adding matches to index..." << std::flush;
    }
    Rcpp::RObject match_table = tar_read(match_table_name);
    Rcpp::IntegerVector seq_idx1 =
      integer_column(match_table, "seq_idx1", match_table_name.get_cstring());
    Rcpp::IntegerVector seq_idx2 =
      integer_column(match_table, "seq_idx2", match_table_name.get_cstring());
    Rcpp::IntegerVector nread1 =
      integer_column(match_table, "nread1", match_table_name.get_cstring());
    Rcpp::IntegerVector nread2 =
      integer_column(match_table, "nread2", match_table_name.get_cstring());
    Rcpp::NumericVector dist =
      numeric_column(match_table, "dist", match_table_name.get_cstring());

    for (R_xlen_t j = 0; j < seq_idx1.size(); ++j) {
      if (Rcpp::IntegerVector::is_na(seq_idx1[j])) continue;
      if (Rcpp::IntegerVector::is_na(seq_idx2[j])) continue;
      if (Rcpp::IntegerVector::is_na(nread1[j])) continue;
      if (Rcpp::IntegerVector::is_na(nread2[j])) continue;
      if (Rcpp::NumericVector::is_na(dist[j])) continue;
      if (Rcpp::traits::is_nan<REALSXP>(dist[j])) continue;
      if (dist[j] > max_dist) continue;
      mid.add_match(
        fwd_map.at(seq_idx1[j]),
        fwd_map.at(seq_idx2[j]),
        nread1[j],
        nread2[j]
      );
    }
    if (verbose) Rcpp::Rcerr << "done." << std::endl;
  }
  if (verbose) Rcpp::Rcerr << "Collecting garbage..." << std::flush;
  R_gc();
  if (verbose) Rcpp::Rcerr << "done." << std::endl;

  // initialize seq_idx_out with the mapped indices
  Rcpp::IntegerVector seq_idx_out(n_seq_idx);
  // nonempty_occurrences indexed by mapped sequence indices
  std::vector<int> mapped_total_occurrences(n_seq_idx);

  for (int i = 0; i < rev_map.size(); ++i) {
    seq_idx_out[i] = i;
    mapped_total_occurrences[i] = total_occurrences[rev_map[i]];
  }

  Rcpp::DataFrame lulu_map = lulu_map_impl(
     seq_idx_out,
     mid,
     mapped_total_occurrences,
     min_abundance_ratio,
     min_cooccurrence_ratio,
     use_mean_abundance_ratio,
     verbose
  );

  // now apply the reverse map
  // seq_idx_out was added directly to the data frame without modification,
  // so our original handle to it is still valid!
  Rcpp::IntegerVector lulu_idx = lulu_map["lulu_idx"];
  for (int i = 0; i < seq_idx_out.size(); ++i) {
    seq_idx_out[i] = rev_map.at(seq_idx_out[i]);
    lulu_idx[i] = rev_map.at(lulu_idx[i]);
  }

  return lulu_map;
}

