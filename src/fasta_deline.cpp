#include "fastq.h"

//' Convert a multiline FASTA file into a single-line FASTA file
//'
//' @param infile (`character(1)`) FASTA file, can be compressed
//' @param outfile (`character(1)`) name of output FASTA file
//' @param compress (`logical(1)`) compress output file
//' @return `character(1)` output file name
//' @export
// [[Rcpp::export]]
std::string fasta_deline(
  const std::string infile,
  const std::string outfile,
  const bool compress = false
) {
  filter_stream_in in_stream;
  open_fastx_in(in_stream, infile);
  if (!in_stream) {
    Rcpp::stop("Error opening FASTA file %s for reading", infile);
  }

  filter_stream_out out_stream;
  open_fastx_out(out_stream, outfile, compress);
  if (!out_stream) {
    Rcpp::stop("Error opening FASTA file %s for writing", outfile);
  }

  std::string header, seq;
  std::getline(in_stream, header);
  while (in_stream) {
    if (header.empty()) {
      continue;
    }
    out_stream << header << '\n';
    while (std::getline(in_stream, seq) && seq[0] != '>') {
      const auto end = seq.find_last_not_of("\r");
      out_stream << seq.substr(0, end + 1);
    }
    out_stream << '\n';
  }
  return outfile;
}
