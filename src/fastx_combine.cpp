#include "fastq.h"

//' Interleave multiple FASTQ files
//' @param infiles (`character()`) input FASTQ files, can be compressed
//' @param outfile (`character(1)`) output FASTQ file
//' @param compress (`logical(1)`) compress output file
//' @return `character(1)` output file name
//' @export
//' @keywords internal
// [[Rcpp::export]]
std::string fastq_combine(
  const std::vector<std::string> & infiles,
  const std::string outfile,
  const bool compress = false
) {
  std::vector<filter_stream_in> in_streams(infiles.size());
  for (size_t i = 0; i < infiles.size(); ++i) {
    open_fastx_in(in_streams[i], infiles[i]);
  }

  filter_stream_out out_stream;
  open_fastx_out(out_stream, outfile, compress);

  std::string header, seq, header2, qual;
  int i = 0;
  bool all_eof = false;
  while (!all_eof) {
    all_eof = true;
    for (size_t j = 0; j < in_streams.size(); ++j) {
      if (in_streams[j].eof()) {
        continue;
      }
      in_streams[j] >> header;
      if (in_streams[j]) {
        all_eof = false;
        if (in_streams[j] >> seq >> header2 >> qual) {
          out_stream << header << '\n'
                     << seq << '\n'
                     << header2 << '\n'
                     << qual << '\n';
        } else {
          Rcpp::stop("Error reading FASTQ file %s", infiles[j]);
        }
      }
    }
  }
  return outfile;
}

//' Interleave multiple FASTA files
//' @param infiles (`character()`) input FASTA files, can be compressed
//' @param outfile (`character(1)`) output FASTA file
//' @param compress (`logical(1)`) compress output file
//' @return `character(1)` output file name
//' @export
//' @keywords internal
// [[Rcpp::export]]
std::string fasta_combine(
  const std::vector<std::string> & infiles,
  const std::string outfile,
  const bool compress = false
) {
  std::vector<filter_stream_in> in_streams(infiles.size());
  for (size_t i = 0; i < infiles.size(); ++i) {
    open_fastx_in(in_streams[i], infiles[i]);
  }

  filter_stream_out out_stream;
  open_fastx_out(out_stream, outfile, compress);

  std::vector<std::string> line(in_streams.size());
  for (size_t i = 0; i < in_streams.size(); ++i) {
    if (!(in_streams[i] >> line[i])) {
      Rcpp::stop("Error reading FASTA file %s", infiles[i]);
    }
    if (line[i][0] != '>') {
      Rcpp::stop("Expected header in FASTA file %s", infiles[i]);
    }
  }
  bool all_eof = false;
  while (!all_eof) {
    all_eof = true;
    for (size_t j = 0; j < in_streams.size(); ++j) {
      if (in_streams[j].eof()) {
        continue;
      }
      all_eof = false;
      out_stream << line[j] << '\n';
      in_streams[j] >> line[j];
      if (line[j][0] == '>') {
        Rcpp::stop("Expected sequence in FASTA file %s", infiles[j]);
      }
      do {
        out_stream << line[j] << '\n';
      } while (in_streams[j] >> line[j] && line[j][0] != '>');
    }
  }
  return outfile;
}
