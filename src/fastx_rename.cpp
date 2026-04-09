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
  std::string line;
  while (std::getline(instream, line)) {
    if (line.empty()) continue;
    if (line[0] != '@') {
      Rcpp::stop("Malformed FASTQ: expected '@' header line");
    }
    if (i >= static_cast<std::size_t>(names.size())) {
      Rcpp::stop("Not enough replacement names for FASTQ records");
    }

    outstream << "@" << names[i++] << "\n";

    std::string seq;
    std::string plus;
    std::string qual;
    if (!std::getline(instream, seq)
        || !std::getline(instream, plus)
        || !std::getline(instream, qual)) {
      Rcpp::stop("Malformed FASTQ: incomplete record");
    }

    outstream << seq << "\n";
    outstream << plus << "\n";
    outstream << qual << "\n";
  }

  if (i != static_cast<std::size_t>(names.size())) {
    Rcpp::stop("Number of replacement names does not match FASTQ records");
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
  std::string line;
  while (std::getline(instream, line)) {
    if (line.empty()) continue;
    if (line[0] == '>') {
      if (i >= static_cast<std::size_t>(names.size())) {
        Rcpp::stop("Not enough replacement names for FASTA records");
      }
      outstream << ">" << names[i++] << "\n";
    } else {
      outstream << line << "\n";
    }
  }

  if (i != static_cast<std::size_t>(names.size())) {
    Rcpp::stop("Number of replacement names does not match FASTA records");
  }

  return outfile;
}
