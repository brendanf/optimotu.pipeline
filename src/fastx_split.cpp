#include "fastq.h"

//' Split a FASTQ file into multiple files
//' @param infile Input FASTQ file
//' @param outfiles Output FASTQ files
//' @param compress Whether to compress the output files
//' @return Output file names
//' @export
//' @keywords internal
// [[Rcpp::export]]
Rcpp::CharacterVector fastq_split(
    std::string infile,
    std::vector<std::string> outfiles,
    bool compress = false) {

  filter_stream_in instream;
  open_fastx_in(instream, infile);

  std::vector<filter_stream_out> outstreams(outfiles.size());
  for (size_t i = 0; i < outfiles.size(); ++i) {
    open_fastx_out(outstreams[i], outfiles[i], compress);
  }

  size_t i = 0;
  std::string header, seq, header2, qual;
  while (instream >> header) {
    if (!(instream >> seq >> header2 >> qual)) {
      Rcpp::stop("Unexpected end of file");
    }
    outstreams[i] << header << "\n"
                  << seq << "\n"
                  << header2 << "\n"
                  << qual << "\n";
    i = (i + 1) % outfiles.size();
  }

  return Rcpp::wrap(outfiles);
}

//' Split a FASTA file into multiple files
//' @param infile Input FASTA file
//' @param outfiles Output FASTA files
//' @param compress Whether to compress the output files
//' @return Output file names
//' @export
//' @keywords internal
// [[Rcpp::export]]
Rcpp::CharacterVector fasta_split(
    std::string infile,
    std::vector<std::string> outfiles,
    bool compress = false) {
  filter_stream_in instream;
  open_fastx_in(instream, infile);

  std::vector<filter_stream_out> outstreams(outfiles.size());
  for (size_t i = 0; i < outfiles.size(); ++i) {
    open_fastx_out(outstreams[i], outfiles[i], compress);
  }

  size_t i = 0;
  if (!instream.good()) {
    Rcpp::stop("Empty file");
  }
  std::string line;
  instream >> line;
  while (instream.good()) {
    if (!boost::starts_with(line, ">")) {
      Rcpp::stop("Expected header line");
    }
    outstreams[i] << line << "\n";
    if (!(instream >> line)) {
      Rcpp::stop("Unexpected end of file");
    }
    if (boost::starts_with(line, ">")) {
      Rcpp::stop("Expected sequence line");
    }
    do {
      outstreams[i] << line << "\n";
    } while (instream >> line && !boost::starts_with(line, ">"));
    i = (i + 1) % outfiles.size();
  }
  return Rcpp::wrap(outfiles);
}
