#include <Rcpp.h>
#include <cstdint>

//' Calculate empirical transition matrix from known true sequences
//'
//' @param true_seq (`character`) Vector of true sequences
//' @param obs_seq (`character`) Vector of observed reads
//' @param obs_qual (`character`) Vector of quality scores for observed reads
//' @param match_idx (`integer`) Indices of the true sequences that match each
//' observed read
//' @param cigar (`character`) CIGAR strings for alignment between each
//' observed read and its corresponding true sequence
//' @param scores (`integer`) set of unique quality scores found in the
//' observed reads
//' @param indel (`logical`) Whether to include indels in the transition
//' matrix
//'
//' Calculates an empirical translation matrix based on pairwise alignments
//' between true sequences (e.g., from a positive control or mock community)
//' and observed reads. In the case of substitutions (i.e., matches and
//' mismatches), there are 16 possible transitions (4 bases each for the true
//' base and the observed base), each of which is represented as one row in the
//' resulting matrix. The rownames are of the form "X2Y", representing
//' "observed cases where X in the true sequence was observed as Y in a read",
//' where both X and Y are "A", "C", "G", or "T". The columns represent the
//' quality scores of the base in the reads.
//'
//' If `indel = TRUE`, the matrix will also include transitions for insertions
//' and deletions. Because each of these represents a lack of correspondence
//' between the true sequence and the observed read, there is some ambiguity
//' in exactly how they should be mapped.
//'
//' For insertions, the inserted base(s) have no corresponding element in the
//' true sequence, but may in theory be mapped to either the previous base or
//' the next base in the true sequence. The matrix includes two rows for
//' these, with names of the form "X_2NY", and "_X2YN", meaning "The position
//' after X in the true sequence had an observed insertion of Y" and "The
//' position before X in the true sequence had an observed insertion of Y",
//' respectively. "N" in this case symbolizes "any base in the position of X",
//' because it is possible that the insertion is directly preceded or followed
//' by a substitution. The column for the insertion is the quality score of the
//' inserted base (i.e., "Y").
//'
//' For deletions, it is clear which base in the true sequence was deleted, but
//' the quality score of the deleted base is not available; however it is
//' expected that the quality score of the preceding or following base might be
//' indicative of the possible presence of a deletion. Thus each deletion
//' occurs twice in the matrix, once for each of the two possible quality
//' scores. The row names for these are "XN2_N" and "NX2N_", meaning in both
//' cases "X in the true sequence was deleted in the observed sequence"; the
//' "N" is a placeholder for the base in the observed sequence that supplies
//' the quality score for the deletion.
//'
//' @return An integer matrix with rows representing transitions from one base
//' to another (or to an insertion or deletion when `indel = TRUE`), columns
//' representing quality scores, and values representing the count of each
//' transition.
//'
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix empirical_transition_matrix(
    Rcpp::CharacterVector true_seq,
    Rcpp::CharacterVector obs_seq,
    Rcpp::CharacterVector obs_qual,
    Rcpp::IntegerVector match_idx,
    Rcpp::CharacterVector cigar,
    Rcpp::IntegerVector scores,
    bool indel = false
) {
  // number of possible transitions from each base
  int base_offset = 4;
  if (indel) {
    // 4 possible transitions, 4 insertions before, 4 insertions after, 1 deletion
    base_offset = 14;
  }
  int n_row = 4 * base_offset;
  Rcpp::IntegerMatrix transition_matrix(n_row, scores.size());
  Rcpp::colnames(transition_matrix) = scores;
  if (indel) {
    Rcpp::rownames(transition_matrix) = Rcpp::CharacterVector::create(
      "A2A", "A2C", "A2G", "A2T", "A_2NA", "A_2NC", "A_2NG", "A_2NT", "_A2AN", "_A2CN", "_A2GN", "_A2TN", "AN2_N", "NA2N_",
      "C2A", "C2C", "C2G", "C2T", "C_2NA", "C_2NC", "C_2NG", "C_2NT", "_C2AN", "_C2CN", "_C2GN", "_C2TN", "CN2_N", "NC2N_",
      "G2A", "G2C", "G2G", "G2T", "G_2NA", "G_2NC", "G_2NG", "G_2NT", "_G2AN", "_G2CN", "_G2GN", "_G2TG", "GN2_N", "NG2N_",
      "T2A", "T2C", "T2G", "T2T", "T_2NA", "T_2NC", "T_2NG", "T_2NT", "_T2AN", "_T2CN", "_T2GN", "_T2TT", "TN2_N", "NT2N_"
    );
  } else {
    Rcpp::rownames(transition_matrix) = Rcpp::CharacterVector::create(
      "A2A", "A2C", "A2G", "A2T",
      "C2A", "C2C", "C2G", "C2T",
      "G2A", "G2C", "G2G", "G2T",
      "T2A", "T2C", "T2G", "T2T"
    );
  }

  int max_score = Rcpp::max(scores) + 33;
  std::vector<int> score_map(max_score + 1, -1);
  for (int i = 0; i < scores.size(); i++) {
    score_map[scores[i] + 33] = i;
  }

  for (int i = 0; i < match_idx.size(); i++) {
    int idx = match_idx[i] - 1; // Convert to zero-based index
    Rcpp::String true_seq_i = true_seq[idx];
    const char * true_base = true_seq_i.get_cstring();
    const char * prev_true = nullptr;


    Rcpp::String obs_seq_i = obs_seq[i];
    const char * obs_base = obs_seq_i.get_cstring();
    Rcpp::String obs_qual_i = obs_qual[i];
    const char* qual = obs_qual_i.get_cstring();
    const char* prev_qual = nullptr;

    int row = -1, col = -1;

    if (cigar[i] == "=") {
      // VSEARCH outputs a single "=" for a complete match.

      while(*true_base && *obs_base) {
        std::uint8_t row;
        if (*true_base == 'A') {
          row = 0;
        } else if (*true_base == 'C') {
          row = base_offset;
        } else if (*true_base == 'G') {
          row = 2 * base_offset;
        } else if (*true_base == 'T') {
          row = 3 * base_offset;
        } else {
          Rcpp::stop("Invalid base in true sequence");
        }
        if (*qual < 33 || *qual >= (int)score_map.size()) {
          Rcpp::stop("Invalid quality score");
        }
        int col = score_map[*qual];
        if (col < 0) {
          Rcpp::stop("Invalid quality score");
        }
        if (col > 2) {
          // col == 2 should be "N"; vsearch does not consider this a mismatch
          // but it is not really a match either.
          transition_matrix(row, col)++;
        }
        true_base++;
        obs_base++;
        qual++;
      }
    } else {
      Rcpp::String cigar_i = cigar[i];
      const char* cigar_op = cigar_i.get_cstring();
      while (*cigar_op && *true_base && *obs_base && *qual) {
        // otherwise expect zero or more digits followed by a letter.
        std::size_t count = 0;
        while (*cigar_op && isdigit(*cigar_op)) {
          count = count * 10 + (*cigar_op - '0');
          cigar_op++;
        }
        // if no digits, default to 1
        if (count == 0) {
          count = 1;
        }
        // check for early end
        if (!*cigar_op) {
          Rcpp::stop("Invalid CIGAR string ends with number: " + cigar[i]);
        }

        for (std::size_t k = 0; k < count; k++) {
          if (*cigar_op == 'M' || *cigar_op == 'X' || *cigar_op == '=') {

            if (*true_base == 'A') {
              row = 0;
            } else if (*true_base == 'C') {
              row = base_offset;
            } else if (*true_base == 'G') {
              row = 2 * base_offset;
            } else if (*true_base == 'T') {
              row = 3 * base_offset;
            } else {
              Rcpp::stop("Invalid base in true sequence");
            }
            if (*obs_base == 'A') {
              row += 0;
            } else if (*obs_base == 'C') {
              row += 1;
            } else if (*obs_base == 'G') {
              row += 2;
            } else if (*obs_base == 'T') {
              row += 3;
            } else if (*obs_base == 'N') {
              true_base++;
              obs_base++;
              qual++;
              continue; // ignore N
            } else if (*obs_base == 0) {
              Rcpp::stop("Unexpended end of string in observed sequence " + std::to_string(i + 1));
            } else {
              Rcpp::stop("Invalid base in observed sequence " + std::to_string(i + 1) + ": " + std::string(1, *obs_base));
            }
            if (*qual < 33 || *qual >= (int)score_map.size()) {
              Rcpp::stop("Invalid quality score");
            }
            col = score_map[*qual];
            if (col < 0 || col >= scores.size() || row < 0 || row >= n_row) {
              Rcpp::stop("Invalid row or column index: row=" +
                std::to_string(row) + ", col=" + std::to_string(col));
            }
            // advance both true and observed sequence
            prev_true = true_base;
            true_base++;
            obs_base++;
            prev_qual = qual;
            qual++;
            transition_matrix(row, col)++;
          } else if (*cigar_op == 'I') {
            if (indel) {
              if (*qual < 33 || *qual >= (int)score_map.size()) {
                Rcpp::stop("Invalid quality score");
              }
              col = score_map[*qual];

              // insertion gets added twice. First after the previous true base
              if (prev_true) {
                if (*prev_true == 'A') {
                  row = 0;
                } else if (*prev_true == 'C') {
                  row = base_offset;
                } else if (*prev_true == 'G') {
                  row = 2 * base_offset;
                } else if (*prev_true == 'T') {
                  row = 3 * base_offset;
                } else {
                  Rcpp::stop("Invalid base in previous observed sequence");
                }
                if (*obs_base == 'A') {
                  row += 4;
                } else if (*obs_base == 'C') {
                  row += 5;
                } else if (*obs_base == 'G') {
                  row += 6;
                } else if (*obs_base == 'T') {
                  row += 7;
                } else if (*obs_base == 'N') {
                  row = -1; // ignore N
                } else if (*obs_base == 0) {
                  Rcpp::stop("Unexpended end of string in observed sequence " +
                    std::to_string(i + 1));
                } else {
                  Rcpp::stop("Invalid base in observed sequence " +
                    std::to_string(i + 1) + ": " + std::string(1, *obs_base));
                }
                if (row >= 0) {
                  if (col < 0 || col >= scores.size() || row >= n_row) {
                    Rcpp::stop("Invalid row or column index: row=" +
                      std::to_string(row) + ", col=" + std::to_string(col));
                  }
                  transition_matrix(row, col)++;
                }
              }
              // now before the next true base
              if (*true_base) {
                if (*true_base == 'A') {
                  row = 0;
                } else if (*true_base == 'C') {
                  row = base_offset;
                } else if (*true_base == 'G') {
                  row = 2 * base_offset;
                } else if (*true_base == 'T') {
                  row = 3 * base_offset;
                } else {
                  Rcpp::stop("Invalid base in true sequence");
                }
                if (*obs_base == 'A') {
                  row += 8;
                } else if (*obs_base == 'C') {
                  row += 9;
                } else if (*obs_base == 'G') {
                  row += 10;
                } else if (*obs_base == 'T') {
                  row += 11;
                } else if (*obs_base == 'N') {
                  row = -1; // ignore N
                } else if (*obs_base == 0) {
                  Rcpp::stop("Unexpended end of string in observed sequence " +
                    std::to_string(i + 1));
                } else {
                  Rcpp::stop("Invalid base in observed sequence " +
                    std::to_string(i + 1) + ": " + std::string(1, *obs_base));
                }
                if (row >= 0) {
                  if (col < 0 || col >= scores.size() || row >= n_row) {
                    Rcpp::stop("Invalid row or column index: row=" +
                      std::to_string(row) + ", col=" + std::to_string(col));
                  }
                  transition_matrix(row, col)++;
                }
              }
            }
            // advanced observed sequence only
            obs_base++;
            prev_qual = qual;
            qual++;
          } else if (*cigar_op == 'D') {
            if (indel) {
              if (*true_base == 'A') {
                row = 0;
              } else if (*true_base == 'C') {
                row = base_offset;
              } else if (*true_base == 'G') {
                row = 2 * base_offset;
              } else if (*true_base == 'T') {
                row = 3 * base_offset;
              } else {
                Rcpp::stop("Invalid base in true sequence");
              }
              row += 12; // first deletion row offset
              if (prev_qual) {
                if (*prev_qual < 33 || *prev_qual >= (int)score_map.size()) {
                  Rcpp::stop("Invalid quality score");
                }
                col = score_map[*prev_qual];
                if (col < 0 || col >= scores.size() || row < 0 || row >= n_row) {
                  Rcpp::stop("Invalid row or column index: row=" +
                    std::to_string(row) + ", col=" + std::to_string(col));
                }
                transition_matrix(row, col)++;
              }
              row++; // second deletion row offset
              if (*qual) {
                if (*qual < 33 || *qual >= (int)score_map.size()) {
                  Rcpp::stop("Invalid quality score");
                }
                col = score_map[*qual];
                if (col < 0 || col >= scores.size() || row < 0 || row >= n_row) {
                  Rcpp::stop("Invalid row or column index: row=" +
                    std::to_string(row) + ", col=" + std::to_string(col));
                }
                transition_matrix(row, col)++;
              }
            }
            // advance true sequence only
            prev_true = true_base;
            true_base++;
          } else {
            Rcpp::stop("Invalid CIGAR operation: " + std::string(1, *cigar_op));
          }
        }
        cigar_op++;
      }
    }
  }
  return transition_matrix;
}
