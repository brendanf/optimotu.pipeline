#include "fastq.h"
#include <cstdint>
#include <array>

//' Get the unique quality scores in a set of FASTQ files
//' @param fastq (`character`) FASTQ file paths, optionally gzipped.
//' @param max_n (`integer`) Maximum number of reads to process.
//' @param offset (`integer`) ASCII offset for quality scores.
//' @return (`integer`) Unique quality scores.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector fastq_qual_bins(
    std::vector<std::string> fastq,
    int max_n = 1e9,
    int offset = 33
) {
  int i = 0;
  std::array<bool, 256> hits;
  hits.fill(false);

  for (const auto & f : fastq) {
    filter_stream_in instream;
    open_fastx_in(instream, f);
    while (!instream.eof()) {
      std::uint8_t c = 0;
      instream.ignore(NO_LIMIT, '\n');
      instream.ignore(NO_LIMIT, '\n');
      instream.ignore(NO_LIMIT, '\n');
      while ((c = instream.get()) != '\n') {
        if (instream.eof()) break;
        if (c == 13) continue; // ignore carriage return (for Windows)
        hits[c] = true;
      }
      ++i;
      if (i == max_n) break;
      Rcpp::checkUserInterrupt();
    }
    if (i == max_n) break;
  }
  Rcpp::IntegerVector out;
  for (int c = 0; c <=255; ++c) {
    if (hits[c]) out.push_back((int)c - offset);
  }
  return out;
}
