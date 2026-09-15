// implementation of LULU algorithm for a long OTU table
// Written without reference to the LULU code, based only on algorithm
// description in the manuscript and documentation.

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include "accessor.h"

// struct to store match information about a pair of sequences
// abund_ratio is either the sum of abundance ratios in all samples where the
//  sequences co-occur (if use_mean_abundance = true) or the minimum abundance
//  ratio (if use_mean_abundance_ratio = false).
// nboth is the number of samples in which they co-occur.
// It is assumed that the first sequence is the potential child, i.e. the one
// which is less prevalent, or if tied then less abundant.
struct match_info
{
  double abund_ratio = 0;
  int nboth = 0;

  match_info(int nread1, int nread2) : abund_ratio(double(nread2) / double(nread1)),
                                       nboth(1) {}

  match_info() {}

  void add_match(int nread1, int nread2, bool use_mean_abundance_ratio)
  {
    double new_abund_ratio = double(nread2) / double(nread1);
    if (nboth == 0)
    {
      abund_ratio = new_abund_ratio;
    }
    else if (use_mean_abundance_ratio)
    {
      abund_ratio += new_abund_ratio;
    }
    else
    {
      abund_ratio = std::min(abund_ratio, new_abund_ratio);
    }
    ++nboth;
  }
};

struct match_info_data
{
  std::map<std::pair<int, int>, match_info> data;
  const bool use_mean;

  match_info_data(bool use_mean_abundance_ratio) : use_mean(use_mean_abundance_ratio) {}

  void add_match(int id1, int id2, int nread1, int nread2)
  {
    // ensure that id1 is the larger ID (i.e., the potential child)
    // this means that we only need to check one direction for matches
    if (id2 > id1)
    {
      std::swap(id1, id2);
      std::swap(nread1, nread2);
    }
    data[std::make_pair(id1, id2)].add_match(nread1, nread2, use_mean);
  }
};

// overload for core implementation
Rcpp::DataFrame lulu_map_impl(
    Rcpp::IntegerVector seq_idx_out,
    match_info_data &match_info,
    std::vector<int> &total_occurrence,
    double min_abundance_ratio = 1.0,
    double min_cooccurrence_ratio = 0.95,
    bool use_mean_abundance_ratio = false,
    int verbose = 0)
{

  std::vector<int> lulu_map(total_occurrence.size(), NA_INTEGER);
  for (int i : seq_idx_out)
  {
    lulu_map[i] = i;
  }

  for (const auto &mi : match_info.data)
  {
    if (verbose > 0)
    {
      Rcpp::Rcerr << "Considering match pair (" << mi.first.first
                  << ", " << mi.first.second << ") with "
                  << total_occurrence[mi.first.first] << " and "
                  << total_occurrence[mi.first.second] << " occurrences and "
                  << mi.second.nboth << " co-occurrences"
                  << std::endl;
    }
    // if the potential child has already been denoised, skip
    if (lulu_map[mi.first.first] != mi.first.first)
    {
      if (verbose > 0)
      {
        Rcpp::Rcerr << "seq " << mi.first.first
                    << " already mapped to seq " << lulu_map[mi.first.first]
                    << "; skipping" << std::endl;
      }
      continue;
    }

    // check the co-occurrence ratio
    if (mi.second.nboth < min_cooccurrence_ratio * total_occurrence[mi.first.first])
    {
      if (verbose > 0)
      {
        Rcpp::Rcerr << "co-occurrence ratio " << mi.second.nboth
                    << " / " << total_occurrence[mi.first.first]
                    << " = " << double(mi.second.nboth) / double(total_occurrence[mi.first.first])
                    << " is less than minimum " << min_cooccurrence_ratio
                    << "; skipping" << std::endl;
      }
      continue;
    }
    else if (verbose > 1)
    {
      Rcpp::Rcerr << "co-occurrence ratio " << mi.second.nboth
                  << " / " << total_occurrence[mi.first.first]
                  << " = " << double(mi.second.nboth) / double(total_occurrence[mi.first.first])
                  << " is greater than or equal to minimum " << min_cooccurrence_ratio
                  << std::endl;
    }

    // check the abundance ratio
    double abundance_ratio = mi.second.abund_ratio;
    if (use_mean_abundance_ratio)
    {
      // the object has accumulated the sum, so we need to divide.
      abundance_ratio /= mi.second.nboth;
    }
    if (abundance_ratio > min_abundance_ratio)
    {
      if (verbose > 1)
      {
        Rcpp::Rcerr << (use_mean_abundance_ratio ? "mean abundance ratio " : "abundance ratio ")
                    << abundance_ratio << " greater than minimum "
                    << min_abundance_ratio << std::endl;
      }
      if (verbose > 0)
      {
        Rcpp::Rcerr << "Mapping child " << mi.first.first
                    << " to parent " << mi.first.second
                    << std::endl;
      }
      lulu_map[mi.first.first] = mi.first.second;
    }
    else
    {
      if (verbose > 0)
      {
        Rcpp::Rcerr << (use_mean_abundance_ratio ? "mean abundance ratio " : "abundance ratio ")
                    << abundance_ratio << " less than or equal to minimum " << min_abundance_ratio
                    << "; skipping" << std::endl;
      }
    }
  }

  Rcpp::IntegerVector lulu_idx_out(seq_idx_out.size());

  for (R_xlen_t i = 0; i < seq_idx_out.size(); i++)
  {
    int j = seq_idx_out[i];
    while (lulu_map[j] != j)
    {
      j = lulu_map[j];
    }
    lulu_idx_out[i] = j;
  }

  return Rcpp::DataFrame::create(
      Rcpp::Named("seq_idx") = seq_idx_out,
      Rcpp::Named("lulu_idx") = lulu_idx_out);
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
    int verbose = 0)
{
  // check that all match_* vectors are the same length
  if (match_id1.size() != match_id2.size() ||
      match_id1.size() != match_nread1.size() ||
      match_id1.size() != match_nread2.size() ||
      match_id1.size() != match_dist.size())
  {
    Rcpp::stop("All match_* vectors must be the same length");
  }
  // check that seq_idx and seq_nsample are the same length
  if (seq_idx.size() != nread.size())
  {
    Rcpp::stop("seq_idx and nread must be the same length");
  }

  // check that the sequences are correctly sorted by occurrence and abundance
  std::vector<int> total_occurrence;
  std::vector<int> total_abundance;
  int n_seq_idx = 0;
  for (int i = 0; i < seq_idx.size(); i++)
  {
    if (seq_idx[i] >= (int)total_occurrence.size())
    {
      total_occurrence.resize(seq_idx[i] + 1, 0);
      total_abundance.resize(seq_idx[i] + 1, 0);
    }
    if (total_occurrence[seq_idx[i]] == 0)
    {
      n_seq_idx++;
    }
    total_occurrence[seq_idx[i]]++;
    total_abundance[seq_idx[i]] += nread[i];
  }

  Rcpp::IntegerVector seq_idx_out(n_seq_idx);
  Rcpp::IntegerVector lulu_idx_out(n_seq_idx);
  int j = 0;

  for (std::size_t i = 0; i < total_occurrence.size(); i++)
  {
    if (total_occurrence[i] > 0)
    {
      seq_idx_out[j] = i;
      j++;
    }
  }

  for (R_xlen_t i = 1; i < seq_idx_out.size(); i++)
  {
    if (total_occurrence[seq_idx_out[i]] > total_occurrence[seq_idx_out[i - 1]])
    {
      Rcpp::stop(
          "seq_idx %d has %d occurences, greater than seq_idx %d with %d.",
          seq_idx_out[i],
          total_occurrence[seq_idx_out[i]],
          seq_idx_out[i - 1],
          total_occurrence[seq_idx_out[i - 1]]);
    }
    if (total_occurrence[seq_idx_out[i]] == total_occurrence[seq_idx_out[i - 1]] &&
        total_abundance[seq_idx_out[i]] > total_abundance[seq_idx_out[i - 1]])
    {
      Rcpp::stop("seq_idx must be sorted by decreasing occurrence, with ties"
                 "broken by decreasing abundance");
    }
  }

  match_info_data match_info(use_mean_abundance_ratio);
  for (int i = 0; i < match_id1.size(); i++)
  {
    if (match_dist[i] <= max_dist)
    {
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
      verbose);
}

Rcpp::RObject tar_read(Rcpp::String name)
{
  // Crew workers leave tar_runtime$meta unset, so tar_read_raw() is
  // forbidden there. read_runtime_target() uses in-memory meta or the
  // current target's subpipeline instead.
  Rcpp::Environment pkg =
      Rcpp::Environment::namespace_env("optimotu.pipeline");
  Rcpp::Function read_runtime_target = pkg["read_runtime_target"];
  return read_runtime_target(name);
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

)
{
  // First read the OTU table(s) to get totals
  std::vector<int> total_occurrences;
  std::vector<size_t> total_abundance;

  for (R_xlen_t i = 0; i < otu_table_names.size(); ++i)
  {
    Rcpp::String otu_table_name(otu_table_names[i]);
    if (verbose)
    {
      Rcpp::Rcerr << "Reading OTU table " << otu_table_name.get_cstring()
                  << "\n  Collecting garbage..." << std::flush;
    }
    R_gc();
    if (verbose)
    {
      Rcpp::Rcerr << "done.\n  Counting occurrences..." << std::flush;
    }
    Rcpp::RObject otu_table = tar_read(otu_table_name);
    Rcpp::IntegerVector seq_idx =
        integer_column(otu_table, "seq_idx", otu_table_name.get_cstring());
    Rcpp::IntegerVector nread =
        integer_column(otu_table, "nread", otu_table_name.get_cstring());

    for (R_xlen_t j = 0; j < seq_idx.size(); ++j)
    {
      if (seq_idx[j] >= (R_xlen_t)total_occurrences.size())
      {
        total_occurrences.resize(seq_idx[j] + 1, 0);
        total_abundance.resize(seq_idx[j] + 1, 0);
      }
      total_occurrences[seq_idx[j]]++;
      total_abundance[seq_idx[j]] += nread[j];
    }
    if (verbose)
      Rcpp::Rcerr << "done." << std::endl;
  }
  if (verbose)
    Rcpp::Rcerr << "Collecting garbage..." << std::flush;
  R_gc();
  if (verbose)
  {
    Rcpp::Rcerr << "done.\nCounting used seq_idx values..." << std::flush;
  }
  // Count the seq indices that actually occur.
  std::size_t n_seq_idx = 0;
  for (int n : total_occurrences)
  {
    if (n > 0)
      ++n_seq_idx;
  }
  if (verbose)
  {
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
  for (int n : total_occurrences)
  {
    if (n > 0)
    {
      rev_map[i] = j;
      nonempty_occurrences[i] = n;
      nonempty_abundance[i] = total_abundance[j];
      ++i;
    }
    ++j;
  }
  if (verbose)
  {
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
      TRUE);

  rev_map = rev_map[order];

  if (verbose)
  {
    Rcpp::Rcerr << "done.\nInverting to form forward map..." << std::flush;
  }

  // invert the rev_map to get the fwd_map
  i = 0;
  for (int seq_idx : rev_map)
  {
    fwd_map[seq_idx] = i++;
  }

  if (verbose)
    Rcpp::Rcerr << "done." << std::endl;

  // Now read the match table(s) to count co-occurrences and relative abundances
  match_info_data mid(use_mean_abundance_ratio);

  for (R_xlen_t i = 0; i < match_table_names.size(); ++i)
  {
    Rcpp::String match_table_name(match_table_names[i]);
    if (verbose)
    {
      Rcpp::Rcerr << "Reading match table " << match_table_name.get_cstring()
                  << "\n  Collecting garbage.." << std::flush;
    }
    R_gc();
    if (verbose)
    {
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

    for (R_xlen_t j = 0; j < seq_idx1.size(); ++j)
    {
      if (Rcpp::IntegerVector::is_na(seq_idx1[j]))
        continue;
      if (Rcpp::IntegerVector::is_na(seq_idx2[j]))
        continue;
      if (Rcpp::IntegerVector::is_na(nread1[j]))
        continue;
      if (Rcpp::IntegerVector::is_na(nread2[j]))
        continue;
      if (Rcpp::NumericVector::is_na(dist[j]))
        continue;
      if (Rcpp::traits::is_nan<REALSXP>(dist[j]))
        continue;
      if (dist[j] > max_dist)
        continue;
      mid.add_match(
          fwd_map.at(seq_idx1[j]),
          fwd_map.at(seq_idx2[j]),
          nread1[j],
          nread2[j]);
    }
    if (verbose)
      Rcpp::Rcerr << "done." << std::endl;
  }
  if (verbose)
    Rcpp::Rcerr << "Collecting garbage..." << std::flush;
  R_gc();
  if (verbose)
    Rcpp::Rcerr << "done." << std::endl;

  // initialize seq_idx_out with the mapped indices
  Rcpp::IntegerVector seq_idx_out(n_seq_idx);
  // nonempty_occurrences indexed by mapped sequence indices
  std::vector<int> mapped_total_occurrences(n_seq_idx);

  for (int i = 0; i < rev_map.size(); ++i)
  {
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
      verbose);

  // now apply the reverse map
  // seq_idx_out was added directly to the data frame without modification,
  // so our original handle to it is still valid!
  Rcpp::IntegerVector lulu_idx = lulu_map["lulu_idx"];
  for (int i = 0; i < seq_idx_out.size(); ++i)
  {
    seq_idx_out[i] = rev_map.at(seq_idx_out[i]);
    lulu_idx[i] = rev_map.at(lulu_idx[i]);
  }

  return lulu_map;
}

// Grain levels for three-level LULU peeling. Narrower grains are smaller.
static constexpr int LULU_GRAIN_BATCH = 0;
static constexpr int LULU_GRAIN_SEQRUN = 1;
static constexpr int LULU_GRAIN_GLOBAL = 2;

static int parse_lulu_scope(Rcpp::String scope)
{
  std::string s(scope.get_cstring());
  if (s == "batch")
    return LULU_GRAIN_BATCH;
  if (s == "seqrun")
    return LULU_GRAIN_SEQRUN;
  if (s == "global")
    return LULU_GRAIN_GLOBAL;
  Rcpp::stop(
      "Unknown LULU scope '%s'; expected 'batch', 'seqrun', or 'global'",
      s.c_str());
  return LULU_GRAIN_GLOBAL;
}

//' Precompute global LULU OTU statistics for scoped mapping
//'
//' Scans OTU table targets one at a time and returns occurrence, abundance,
//' and partition grain for each OTU, sorted in the same parent/child rank
//' order used by [lulu_map_lowmem_impl()].
//'
//' @param otu_table_names (`character`) resolved OTU table target names.
//' @param seqrun_ids (`integer`) parallel to `otu_table_names`; same id means
//'   the tables belong to the same sequencing-run stem.
//' @param verbose (`integer`) verbosity level.
//' @returns a `data.frame` with columns `seq_idx`, `occurrence`, `abundance`,
//'   and `grain` (`0` = batch, `1` = seqrun, `2` = global), sorted most
//'   parent-like first.
//'
// [[Rcpp::export]]
Rcpp::DataFrame lulu_otu_stats_impl(
    Rcpp::CharacterVector otu_table_names,
    Rcpp::IntegerVector seqrun_ids,
    int verbose = 0)
{
  if (otu_table_names.size() != seqrun_ids.size())
  {
    Rcpp::stop("otu_table_names and seqrun_ids must have the same length");
  }

  std::vector<int> total_occurrences;
  std::vector<size_t> total_abundance;
  std::vector<int> first_seqrun;
  std::vector<int> n_tables;
  std::vector<char> multi_seqrun;
  std::vector<int> last_table;

  for (R_xlen_t i = 0; i < otu_table_names.size(); ++i)
  {
    if (Rcpp::IntegerVector::is_na(seqrun_ids[i]))
    {
      Rcpp::stop("seqrun_ids must not contain NA");
    }
    Rcpp::String otu_table_name(otu_table_names[i]);
    if (verbose)
    {
      Rcpp::Rcerr << "Reading OTU table " << otu_table_name.get_cstring()
                  << "\n  Collecting garbage..." << std::flush;
    }
    R_gc();
    if (verbose)
    {
      Rcpp::Rcerr << "done.\n  Counting occurrences..." << std::flush;
    }
    Rcpp::RObject otu_table = tar_read(otu_table_name);
    Rcpp::IntegerVector seq_idx =
        integer_column(otu_table, "seq_idx", otu_table_name.get_cstring());
    Rcpp::IntegerVector nread =
        integer_column(otu_table, "nread", otu_table_name.get_cstring());

    for (R_xlen_t j = 0; j < seq_idx.size(); ++j)
    {
      int s = seq_idx[j];
      if (s < 0)
        continue;
      if (s >= (int)total_occurrences.size())
      {
        std::size_t new_size = (std::size_t)s + 1;
        total_occurrences.resize(new_size, 0);
        total_abundance.resize(new_size, 0);
        first_seqrun.resize(new_size, -1);
        n_tables.resize(new_size, 0);
        multi_seqrun.resize(new_size, 0);
        last_table.resize(new_size, -1);
      }
      total_occurrences[s]++;
      total_abundance[s] += nread[j];
      if (last_table[s] != (int)i)
      {
        last_table[s] = (int)i;
        if (first_seqrun[s] < 0)
        {
          first_seqrun[s] = seqrun_ids[i];
          n_tables[s] = 1;
        }
        else
        {
          n_tables[s]++;
          if (first_seqrun[s] != seqrun_ids[i])
            multi_seqrun[s] = 1;
        }
      }
    }
    if (verbose)
      Rcpp::Rcerr << "done." << std::endl;
  }

  std::size_t n_seq_idx = 0;
  for (int n : total_occurrences)
  {
    if (n > 0)
      ++n_seq_idx;
  }

  Rcpp::IntegerVector rev_map(n_seq_idx);
  Rcpp::IntegerVector nonempty_occurrences(n_seq_idx);
  Rcpp::IntegerVector nonempty_abundance(n_seq_idx);
  Rcpp::IntegerVector nonempty_grain(n_seq_idx);
  Rcpp::IntegerVector order(n_seq_idx);

  int i = 0, j = 0;
  for (int n : total_occurrences)
  {
    if (n > 0)
    {
      rev_map[i] = j;
      nonempty_occurrences[i] = n;
      nonempty_abundance[i] = (int)total_abundance[j];
      if (multi_seqrun[j])
        nonempty_grain[i] = LULU_GRAIN_GLOBAL;
      else if (n_tables[j] > 1)
        nonempty_grain[i] = LULU_GRAIN_SEQRUN;
      else
        nonempty_grain[i] = LULU_GRAIN_BATCH;
      ++i;
    }
    ++j;
  }

  R_orderVector(
      INTEGER(order),
      n_seq_idx,
      Rf_lang2(nonempty_occurrences, nonempty_abundance),
      FALSE,
      TRUE);

  rev_map = rev_map[order];
  nonempty_occurrences = nonempty_occurrences[order];
  nonempty_abundance = nonempty_abundance[order];
  nonempty_grain = nonempty_grain[order];

  return Rcpp::DataFrame::create(
      Rcpp::Named("seq_idx") = rev_map,
      Rcpp::Named("occurrence") = nonempty_occurrences,
      Rcpp::Named("abundance") = nonempty_abundance,
      Rcpp::Named("grain") = nonempty_grain);
}

//' In-memory variant of [lulu_otu_stats_impl()] for tests
//'
//' @param otu_tables (`list`) of OTU `data.frame`s with `seq_idx` and `nread`.
//' @param seqrun_ids (`integer`) parallel to `otu_tables`.
//' @param verbose (`integer`) verbosity level.
//' @returns same structure as [lulu_otu_stats_impl()].
//'
// [[Rcpp::export]]
Rcpp::DataFrame lulu_otu_stats_dfs_impl(
    Rcpp::List otu_tables,
    Rcpp::IntegerVector seqrun_ids,
    int verbose = 0)
{
  if (otu_tables.size() != seqrun_ids.size())
  {
    Rcpp::stop("otu_tables and seqrun_ids must have the same length");
  }

  std::vector<int> total_occurrences;
  std::vector<size_t> total_abundance;
  std::vector<int> first_seqrun;
  std::vector<int> n_tables;
  std::vector<char> multi_seqrun;
  std::vector<int> last_table;

  for (R_xlen_t i = 0; i < otu_tables.size(); ++i)
  {
    if (Rcpp::IntegerVector::is_na(seqrun_ids[i]))
    {
      Rcpp::stop("seqrun_ids must not contain NA");
    }
    Rcpp::RObject otu_table = otu_tables[i];
    std::string name = "otu_tables[[" + std::to_string(i + 1) + "]]";
    Rcpp::IntegerVector seq_idx =
        integer_column(otu_table, "seq_idx", name.c_str());
    Rcpp::IntegerVector nread =
        integer_column(otu_table, "nread", name.c_str());

    for (R_xlen_t j = 0; j < seq_idx.size(); ++j)
    {
      int s = seq_idx[j];
      if (s < 0)
        continue;
      if (s >= (int)total_occurrences.size())
      {
        std::size_t new_size = (std::size_t)s + 1;
        total_occurrences.resize(new_size, 0);
        total_abundance.resize(new_size, 0);
        first_seqrun.resize(new_size, -1);
        n_tables.resize(new_size, 0);
        multi_seqrun.resize(new_size, 0);
        last_table.resize(new_size, -1);
      }
      total_occurrences[s]++;
      total_abundance[s] += nread[j];
      if (last_table[s] != (int)i)
      {
        last_table[s] = (int)i;
        if (first_seqrun[s] < 0)
        {
          first_seqrun[s] = seqrun_ids[i];
          n_tables[s] = 1;
        }
        else
        {
          n_tables[s]++;
          if (first_seqrun[s] != seqrun_ids[i])
            multi_seqrun[s] = 1;
        }
      }
    }
  }

  std::size_t n_seq_idx = 0;
  for (int n : total_occurrences)
  {
    if (n > 0)
      ++n_seq_idx;
  }

  Rcpp::IntegerVector rev_map(n_seq_idx);
  Rcpp::IntegerVector nonempty_occurrences(n_seq_idx);
  Rcpp::IntegerVector nonempty_abundance(n_seq_idx);
  Rcpp::IntegerVector nonempty_grain(n_seq_idx);
  Rcpp::IntegerVector order(n_seq_idx);

  int i = 0, j = 0;
  for (int n : total_occurrences)
  {
    if (n > 0)
    {
      rev_map[i] = j;
      nonempty_occurrences[i] = n;
      nonempty_abundance[i] = (int)total_abundance[j];
      if (multi_seqrun[j])
        nonempty_grain[i] = LULU_GRAIN_GLOBAL;
      else if (n_tables[j] > 1)
        nonempty_grain[i] = LULU_GRAIN_SEQRUN;
      else
        nonempty_grain[i] = LULU_GRAIN_BATCH;
      ++i;
    }
    ++j;
  }

  R_orderVector(
      INTEGER(order),
      n_seq_idx,
      Rf_lang2(nonempty_occurrences, nonempty_abundance),
      FALSE,
      TRUE);

  rev_map = rev_map[order];
  nonempty_occurrences = nonempty_occurrences[order];
  nonempty_abundance = nonempty_abundance[order];
  nonempty_grain = nonempty_grain[order];

  return Rcpp::DataFrame::create(
      Rcpp::Named("seq_idx") = rev_map,
      Rcpp::Named("occurrence") = nonempty_occurrences,
      Rcpp::Named("abundance") = nonempty_abundance,
      Rcpp::Named("grain") = nonempty_grain);
}

// Build dense fwd_map / occurrence / grain / rev_map from stats rows.
// stats rows must be in parent-rank order (most parent-like first).
static void lulu_stats_maps(
    Rcpp::DataFrame stats,
    std::vector<int> &fwd_map,
    std::vector<int> &occurrence,
    std::vector<int> &grain,
    Rcpp::IntegerVector &rev_map,
    std::vector<int> &mapped_occurrence)
{
  Rcpp::IntegerVector seq_idx = stats["seq_idx"];
  Rcpp::IntegerVector occ = stats["occurrence"];
  Rcpp::IntegerVector gr = stats["grain"];
  R_xlen_t n = seq_idx.size();
  rev_map = Rcpp::IntegerVector(n);
  mapped_occurrence.assign(n, 0);

  int max_seq = 0;
  for (R_xlen_t i = 0; i < n; ++i)
  {
    if (seq_idx[i] > max_seq)
      max_seq = seq_idx[i];
  }
  fwd_map.assign((std::size_t)max_seq + 1, NA_INTEGER);
  occurrence.assign((std::size_t)max_seq + 1, 0);
  grain.assign((std::size_t)max_seq + 1, -1);

  for (R_xlen_t i = 0; i < n; ++i)
  {
    int s = seq_idx[i];
    fwd_map[s] = (int)i;
    occurrence[s] = occ[i];
    grain[s] = gr[i];
    rev_map[i] = s;
    mapped_occurrence[i] = occ[i];
  }
}

static void lulu_add_scoped_match(
    int seq1,
    int seq2,
    int nread1,
    int nread2,
    const std::vector<int> &fwd_map,
    const std::vector<int> &occurrence,
    const std::vector<int> &grain,
    int scope,
    bool skip_nested_parent,
    double min_abundance_ratio,
    match_info_data &mid,
    std::vector<int> &singleton_parent)
{
  if (seq1 < 0 || seq2 < 0)
    return;
  if (seq1 >= (int)fwd_map.size() || seq2 >= (int)fwd_map.size())
    return;
  int m1 = fwd_map[seq1];
  int m2 = fwd_map[seq2];
  if (Rcpp::IntegerVector::is_na(m1) || Rcpp::IntegerVector::is_na(m2))
    return;
  if (m1 == m2)
    return;

  int child_m, parent_m, child_n, parent_n, child_seq, parent_seq;
  if (m1 > m2)
  {
    child_m = m1;
    parent_m = m2;
    child_n = nread1;
    parent_n = nread2;
    child_seq = seq1;
    parent_seq = seq2;
  }
  else
  {
    child_m = m2;
    parent_m = m1;
    child_n = nread2;
    parent_n = nread1;
    child_seq = seq2;
    parent_seq = seq1;
  }

  if (grain[child_seq] != scope)
    return;
  if (skip_nested_parent && grain[parent_seq] >= 0 &&
      grain[parent_seq] < grain[child_seq])
    return;

  if (occurrence[child_seq] == 1)
  {
    double abund_ratio = double(parent_n) / double(child_n);
    if (abund_ratio > min_abundance_ratio)
    {
      if (Rcpp::IntegerVector::is_na(singleton_parent[child_m]) ||
          parent_m < singleton_parent[child_m])
      {
        singleton_parent[child_m] = parent_m;
      }
    }
    return;
  }

  mid.add_match(child_m, parent_m, child_n, parent_n);
}

static void lulu_ingest_match_table(
    Rcpp::RObject match_table,
    const char *name,
    double max_dist,
    const std::vector<int> &fwd_map,
    const std::vector<int> &occurrence,
    const std::vector<int> &grain,
    int scope,
    bool skip_nested_parent,
    double min_abundance_ratio,
    match_info_data &mid,
    std::vector<int> &singleton_parent)
{
  Rcpp::IntegerVector seq_idx1 = integer_column(match_table, "seq_idx1", name);
  Rcpp::IntegerVector seq_idx2 = integer_column(match_table, "seq_idx2", name);
  Rcpp::IntegerVector nread1 = integer_column(match_table, "nread1", name);
  Rcpp::IntegerVector nread2 = integer_column(match_table, "nread2", name);
  Rcpp::NumericVector dist = numeric_column(match_table, "dist", name);

  for (R_xlen_t j = 0; j < seq_idx1.size(); ++j)
  {
    if (Rcpp::IntegerVector::is_na(seq_idx1[j]))
      continue;
    if (Rcpp::IntegerVector::is_na(seq_idx2[j]))
      continue;
    if (Rcpp::IntegerVector::is_na(nread1[j]))
      continue;
    if (Rcpp::IntegerVector::is_na(nread2[j]))
      continue;
    if (Rcpp::NumericVector::is_na(dist[j]))
      continue;
    if (Rcpp::traits::is_nan<REALSXP>(dist[j]))
      continue;
    if (dist[j] > max_dist)
      continue;
    lulu_add_scoped_match(
        seq_idx1[j],
        seq_idx2[j],
        nread1[j],
        nread2[j],
        fwd_map,
        occurrence,
        grain,
        scope,
        skip_nested_parent,
        min_abundance_ratio,
        mid,
        singleton_parent);
  }
}

static Rcpp::DataFrame lulu_map_scoped_finish(
    Rcpp::IntegerVector &rev_map,
    std::vector<int> &mapped_occurrence,
    const std::vector<int> &grain_by_seq,
    int scope,
    match_info_data &mid,
    std::vector<int> &singleton_parent,
    double min_abundance_ratio,
    double min_cooccurrence_ratio,
    bool use_mean_abundance_ratio,
    int verbose)
{
  R_xlen_t n_seq_idx = rev_map.size();
  // Immediate parents only (no full path-compression across grains). Within-
  // job chains are compressed below so a batch-restricted child can skip
  // through a batch-restricted intermediate parent.
  std::vector<int> lulu_map(n_seq_idx);
  for (R_xlen_t i = 0; i < n_seq_idx; ++i)
    lulu_map[i] = (int)i;

  for (const auto &mi : mid.data)
  {
    int child = mi.first.first;
    int parent = mi.first.second;
    if (verbose > 0)
    {
      Rcpp::Rcerr << "Considering match pair (" << child << ", " << parent
                  << ") with " << mapped_occurrence[child] << " and "
                  << mapped_occurrence[parent]
                  << " occurrences and " << mi.second.nboth
                  << " co-occurrences" << std::endl;
    }
    if (lulu_map[child] != child)
      continue;
    if (mi.second.nboth <
        min_cooccurrence_ratio * mapped_occurrence[child])
      continue;
    double abundance_ratio = mi.second.abund_ratio;
    if (use_mean_abundance_ratio)
      abundance_ratio /= mi.second.nboth;
    if (abundance_ratio > min_abundance_ratio)
    {
      if (verbose > 0)
      {
        Rcpp::Rcerr << "Mapping child " << child << " to parent " << parent
                    << std::endl;
      }
      lulu_map[child] = parent;
    }
  }

  for (R_xlen_t i = 0; i < n_seq_idx; ++i)
  {
    if (!Rcpp::IntegerVector::is_na(singleton_parent[i]) &&
        lulu_map[i] == (int)i)
    {
      lulu_map[i] = singleton_parent[i];
    }
  }

  std::vector<int> seq_out;
  std::vector<int> lulu_out;
  seq_out.reserve(n_seq_idx);
  lulu_out.reserve(n_seq_idx);
  for (R_xlen_t i = 0; i < n_seq_idx; ++i)
  {
    int child_seq = rev_map[i];
    if (grain_by_seq[child_seq] != scope)
      continue;
    if (lulu_map[i] == (int)i)
      continue;
    int j = lulu_map[i];
    while (lulu_map[j] != j)
      j = lulu_map[j];
    seq_out.push_back(child_seq);
    lulu_out.push_back(rev_map[j]);
  }

  return Rcpp::DataFrame::create(
      Rcpp::Named("seq_idx") = seq_out,
      Rcpp::Named("lulu_idx") = lulu_out);
}

//' Scoped LULU mapping from streamed match-table targets
//'
//' Decides parents only for OTUs whose partition grain matches `scope`.
//' Returns sparse non-identity rows. Requires precomputed stats from
//' [lulu_otu_stats_impl()].
//'
//' @param stats (`data.frame`) output of [lulu_otu_stats_impl()].
//' @param match_table_names (`character`) match table target names to stream.
//' @param scope (`character`) `"batch"`, `"seqrun"`, or `"global"`.
//' @inheritParams lulu_map_impl
//'
// [[Rcpp::export]]
Rcpp::DataFrame lulu_map_scoped_impl(
    Rcpp::DataFrame stats,
    Rcpp::CharacterVector match_table_names,
    Rcpp::String scope,
    double max_dist,
    double min_abundance_ratio = 1.0,
    double min_cooccurrence_ratio = 0.95,
    bool use_mean_abundance_ratio = false,
    int verbose = 0)
{
  int scope_i = parse_lulu_scope(scope);
  bool skip_nested_parent = min_cooccurrence_ratio >= 1.0;

  std::vector<int> fwd_map;
  std::vector<int> occurrence;
  std::vector<int> grain;
  Rcpp::IntegerVector rev_map;
  std::vector<int> mapped_occurrence;
  lulu_stats_maps(
      stats, fwd_map, occurrence, grain, rev_map, mapped_occurrence);

  match_info_data mid(use_mean_abundance_ratio);
  std::vector<int> singleton_parent(rev_map.size(), NA_INTEGER);

  for (R_xlen_t i = 0; i < match_table_names.size(); ++i)
  {
    Rcpp::String match_table_name(match_table_names[i]);
    if (verbose)
    {
      Rcpp::Rcerr << "Reading match table " << match_table_name.get_cstring()
                  << "\n  Collecting garbage.." << std::flush;
    }
    R_gc();
    if (verbose)
    {
      Rcpp::Rcerr << "done.\n  Adding matches to index..." << std::flush;
    }
    Rcpp::RObject match_table = tar_read(match_table_name);
    lulu_ingest_match_table(
        match_table,
        match_table_name.get_cstring(),
        max_dist,
        fwd_map,
        occurrence,
        grain,
        scope_i,
        skip_nested_parent,
        min_abundance_ratio,
        mid,
        singleton_parent);
    if (verbose)
      Rcpp::Rcerr << "done." << std::endl;
  }
  if (verbose)
    Rcpp::Rcerr << "Collecting garbage..." << std::flush;
  R_gc();
  if (verbose)
    Rcpp::Rcerr << "done." << std::endl;

  return lulu_map_scoped_finish(
      rev_map,
      mapped_occurrence,
      grain,
      scope_i,
      mid,
      singleton_parent,
      min_abundance_ratio,
      min_cooccurrence_ratio,
      use_mean_abundance_ratio,
      verbose);
}

//' Scoped LULU mapping from in-memory match tables
//'
//' Test/direct-use variant of [lulu_map_scoped_impl()] that takes a list of
//' match-table `data.frame`s instead of target names.
//'
//' @param stats (`data.frame`) output of [lulu_otu_stats_impl()] or an
//'   equivalent table.
//' @param match_tables (`list`) of match `data.frame`s.
//' @param scope (`character`) `"batch"`, `"seqrun"`, or `"global"`.
//' @inheritParams lulu_map_impl
//'
// [[Rcpp::export]]
Rcpp::DataFrame lulu_map_scoped_dfs_impl(
    Rcpp::DataFrame stats,
    Rcpp::List match_tables,
    Rcpp::String scope,
    double max_dist,
    double min_abundance_ratio = 1.0,
    double min_cooccurrence_ratio = 0.95,
    bool use_mean_abundance_ratio = false,
    int verbose = 0)
{
  int scope_i = parse_lulu_scope(scope);
  bool skip_nested_parent = min_cooccurrence_ratio >= 1.0;

  std::vector<int> fwd_map;
  std::vector<int> occurrence;
  std::vector<int> grain;
  Rcpp::IntegerVector rev_map;
  std::vector<int> mapped_occurrence;
  lulu_stats_maps(
      stats, fwd_map, occurrence, grain, rev_map, mapped_occurrence);

  match_info_data mid(use_mean_abundance_ratio);
  std::vector<int> singleton_parent(rev_map.size(), NA_INTEGER);

  for (R_xlen_t i = 0; i < match_tables.size(); ++i)
  {
    Rcpp::RObject match_table = match_tables[i];
    std::string name = "match_tables[[" + std::to_string(i + 1) + "]]";
    lulu_ingest_match_table(
        match_table,
        name.c_str(),
        max_dist,
        fwd_map,
        occurrence,
        grain,
        scope_i,
        skip_nested_parent,
        min_abundance_ratio,
        mid,
        singleton_parent);
  }

  return lulu_map_scoped_finish(
      rev_map,
      mapped_occurrence,
      grain,
      scope_i,
      mid,
      singleton_parent,
      min_abundance_ratio,
      min_cooccurrence_ratio,
      use_mean_abundance_ratio,
      verbose);
}

//' Combine sparse scoped LULU maps into a full parent map
//'
//' Starts from identity for every OTU in `stats`, overlays non-identity rows
//' from `sparse_maps`, then path-compresses to roots.
//'
//' @param stats (`data.frame`) with column `seq_idx` (from
//'   [lulu_otu_stats_impl()]).
//' @param sparse_maps (`list`) of `data.frame`s with `seq_idx` and `lulu_idx`.
//' @returns a full `data.frame` with `seq_idx` and `lulu_idx`.
//'
// [[Rcpp::export]]
Rcpp::DataFrame lulu_map_combine_impl(
    Rcpp::DataFrame stats,
    Rcpp::List sparse_maps)
{
  Rcpp::IntegerVector seq_idx = stats["seq_idx"];
  R_xlen_t n = seq_idx.size();
  if (n == 0)
  {
    return Rcpp::DataFrame::create(
        Rcpp::Named("seq_idx") = Rcpp::IntegerVector(),
        Rcpp::Named("lulu_idx") = Rcpp::IntegerVector());
  }

  int max_seq = 0;
  for (R_xlen_t i = 0; i < n; ++i)
  {
    if (seq_idx[i] > max_seq)
      max_seq = seq_idx[i];
  }

  std::vector<int> parent((std::size_t)max_seq + 1, NA_INTEGER);
  for (R_xlen_t i = 0; i < n; ++i)
    parent[seq_idx[i]] = seq_idx[i];

  for (R_xlen_t i = 0; i < sparse_maps.size(); ++i)
  {
    Rcpp::DataFrame m = Rcpp::as<Rcpp::DataFrame>(sparse_maps[i]);
    if (m.nrows() == 0)
      continue;
    Rcpp::IntegerVector s = m["seq_idx"];
    Rcpp::IntegerVector l = m["lulu_idx"];
    for (R_xlen_t j = 0; j < s.size(); ++j)
    {
      if (Rcpp::IntegerVector::is_na(s[j]) ||
          Rcpp::IntegerVector::is_na(l[j]))
        continue;
      if (s[j] < 0 || s[j] > max_seq)
        continue;
      if (Rcpp::IntegerVector::is_na(parent[s[j]]))
        continue;
      parent[s[j]] = l[j];
    }
  }

  Rcpp::IntegerVector lulu_idx(n);
  for (R_xlen_t i = 0; i < n; ++i)
  {
    int j = seq_idx[i];
    // Path-compress with cycle guard
    int guard = 0;
    while (parent[j] != j)
    {
      j = parent[j];
      if (++guard > max_seq + 1)
        Rcpp::stop("Cycle detected while combining LULU maps");
    }
    lulu_idx[i] = j;
  }

  return Rcpp::DataFrame::create(
      Rcpp::Named("seq_idx") = seq_idx,
      Rcpp::Named("lulu_idx") = lulu_idx);
}
