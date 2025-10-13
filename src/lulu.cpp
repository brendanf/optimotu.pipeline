// implementation of LULU algorithm for a long OTU table
// Written without reference to the LULU code, based only on algorithm
// description in the manuscript and documentation.

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>


// struct to store match information
// id1 and id2 are the sequence IDs
// nread is a vector of pairs of (nread1, nread2)
// nsample1 and nsample2 are the number of samples in which the two sequences
// each occur
struct match_info {
  int id1, id2;
  std::vector<std::pair<int, int>> nread;
  int nsample1, nsample2;

  match_info(
    int id1,
    int id2,
    int nread1,
    int nread2,
    const std::vector<int> & total_occurrence,
    const std::vector<int> & total_abundance
  ) : id1(id1), id2(id2), nsample1(total_occurrence[id1]),
     nsample2(total_occurrence[id2]) {
      if (nsample1 > nsample2 || (nsample1 == nsample2 && id1 > id2)) {
        std::swap(this->id1, this->id2);
        std::swap(this->nsample1, this->nsample2);
        nread.push_back(std::make_pair(nread2, nread1));
      } else {
        nread.push_back(std::make_pair(nread1, nread2));
      }
  }

  match_info(
    int id1,
    int id2,
    int nread1,
    int nread2,
    int nsample1,
    int nsample2
  ) : id1(id1), id2(id2), nsample1(nsample1), nsample2(nsample2) {
    nread.push_back(std::make_pair(nread1, nread2));
  }

  match_info() : id1(0), id2(0), nsample1(0), nsample2(0) {}

  bool operator<(const match_info& other) const {
    if (nsample1 < other.nsample1) {
      return true;
    } else if (nsample1 > other.nsample1) {
      return false;
    } else if (nsample2 > other.nsample2) {
      return true;
    } else if (nsample2 < other.nsample2) {
      return false;
    } else if (id1 < other.id1) {
      return true;
    } else if (id1 > other.id1) {
      return false;
    } else {
      return id2 < other.id2;
    }
  }
};

struct match_info_data {
  std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> data;


  void add_match(int id1, int id2, int nread1, int nread2) {
    // ensure that id1 is the larger ID (i.e., the potential child)
    // this means that we only need to check one direction for matches
    if (id2 > id1) {
      std::swap(id1, id2);
      std::swap(nread1, nread2);
    }
    data[std::make_pair(id1, id2)].push_back(std::make_pair(nread1, nread2));
  }
};

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
  for (int i = 0; i < total_occurrence.size(); i++) {
    if (total_occurrence[i] > 0) {
      seq_idx_out[j] = i;
      j++;
    }
  }

  for (std::size_t i = 1; i < seq_idx_out.size(); i++) {
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


  match_info_data match_info;
  for (int i = 0; i < match_id1.size(); i++) {
    if (match_dist[i] <= max_dist) {
      match_info.add_match(match_id1[i], match_id2[i], match_nread1[i], match_nread2[i]);
    }
  }

  std::vector<int> lulu_map(total_occurrence.size(), NA_INTEGER);
  for (int i : seq_idx_out) {
    lulu_map[i] = i;
  }

  for (const auto & mi : match_info.data) {
    if (verbose > 0) {
      Rcpp::Rcerr << "Considering match pair (" << mi.first.first
                  << ", " << mi.first.second << ") with "
                  << mi.second.size() << " co-occurrences"
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
    // mi.second is the vector of (nread1, nread2) pairs; its size is the
    // number of co-occurrences
    if (mi.second.size() < min_cooccurrence_ratio * total_occurrence[mi.first.first]) {
      if (verbose > 0) {
        Rcpp::Rcerr << "co-occurrence ratio " << mi.second.size()
                    << " / " << total_occurrence[mi.first.first]
                    << " = " << double(mi.second.size()) / double(total_occurrence[mi.first.first])
                    << " is less than minimum " << min_cooccurrence_ratio
                    << "; skipping" << std::endl;
      }
      continue;
    }

    // check the abundance ratio
    double abundance_ratio;
    if (use_mean_abundance_ratio) {
      abundance_ratio = 0;
      for (const auto & nread : mi.second) {
        abundance_ratio += double(nread.second) / double(nread.first);
      }
      abundance_ratio /= mi.second.size();
    } else {
      for (const auto & nread : mi.second) {
        abundance_ratio = double(nread.second) / double(nread.first);
        if (abundance_ratio <= min_abundance_ratio) {
          break;
        }
      }
    }
    if (abundance_ratio > min_abundance_ratio) {
      if (verbose > 0) {
        Rcpp::Rcerr << "Mapping child " << mi.first.first
                    << " to parent " << mi.first.second
                    << std::endl;
      }
      lulu_map[mi.first.first] = mi.first.second;
    } else {
      if (verbose > 0) {
        Rcpp::Rcerr << (use_mean_abundance_ratio ? "mean abundance ratio " : "abundance ratio ")
                    << abundance_ratio << " less than minimum " << min_abundance_ratio
                    << "; skipping" << std::endl;
      }
    }
  }

  for (int i = 0; i < seq_idx_out.size(); i++) {
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

