#include "fastq.h"


//' @describeIn fastx_rename Rename reads in a FASTQ file
//' @param infile (`character`) FAST(A/Q) file path, optionally gzipped.
//' @param names (`character`) New read names.
//' @param outfile (`character`) FAST(A/Q) file path, optionally gzipped.
//' @return (`character`) Output file path.
//' @export
// [[Rcpp::export]]
std::string fastq_rename(
    std::string infile,
    Rcpp::CharacterVector names,
    std::string outfile
) {
  filter_stream_in instream;
  open_fastx_in(instream, infile);

  filter_stream_out outstream;
  open_fastx_out(outstream, outfile, false);

  std::size_t i = 0;
  while (instream.get()) {
    std::string line;
    instream >> line;
    if (line.size() == 0) break;
    outstream << "@" << names[i++] << "\n";
    outstream << instream.rdbuf() << "\n";
    outstream << instream.rdbuf() << "\n";
    outstream << instream.rdbuf() << "\n";
  }

  return outfile;
}

//' @describeIn fastx_rename Rename reads in a FASTA file
//' @export
// [[Rcpp::export]]
std::string fasta_rename(
    std::string infile,
    Rcpp::CharacterVector names,
    std::string outfile
) {
  filter_stream_in instream;
  open_fastx_in(instream, infile);

  filter_stream_out outstream;
  open_fastx_out(outstream, outfile, false);
  if (!outstream) {
    Rcpp::stop("Failed to open output file");
  }

  std::size_t i = 0;
  while (instream.get()) {
    std::string line;
    instream >> line;
    if (line.size() == 0) break;
    outstream << ">" << names[i++] << "\n";
    outstream << instream.rdbuf() << "\n";
  }

  return outfile;
}
