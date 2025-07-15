#include <Rcpp.h>
#include "fastq.h"


//' Subsample from a FASTQ file
//'
//' This function subsamples an (optionally gzipped) FASTQ file by randomly
//' selecting a specified number of reads. The file is streamed and only
//' `denominator` reads are stored in memory at any time, which allows it to
//' handle very large files without running out of memory. Sampling is
//' performed with reservoir sampling.
//'
//' Additionally, it guarantees that, if the same file is resampled with the
//' same seed, but a larger number of reads, the result from the first sampling
//' will be a subset of the new, larger sampling.
//'
//' When the number of sequences in the input file is not divisible by the
//' denominator, the remainder sequences will be sampled as if there were enough
//' additional sequences to fill the denominator, so that the actual number of
//' sampled sequences may differ from (file size) * `numerator` / `denominator`
//' by up to max(`numerator`, `denominator`/2).
//'
//' @param file (`character` string) the path to the FASTQ file.
//' @param numerator (`integer` scalar) numerator of the fraction of reads to
//' sample.
//' @param denominator (`integer` scalar) denominator of the fraction of reads
//' to sample.
//' @param output (`character` string) the path to the output file. If it ends
//' in ".gz", the output will be gzipped.
//' @param rename (`logical` scalar) whether to rename the reads in the output;
//' if `TRUE`, the read names will be replaced with a hexadecimal sequential
//' number
//'
//' @return the output file name
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector fastq_sample(
    const std::string& file,
    const int numerator,
    const int denominator,
    const std::string& output,
    bool rename = false
) {

  if (numerator <= 0 || denominator <= 0 || numerator > denominator) {
    Rcpp::stop("Numerator and denominator must be positive integers, and"
                 " numerator must be less than or equal to denominator.");
  }
  if (file.empty()) {
    Rcpp::stop("Input file name cannot be empty.");
  }
  if (output.empty()) {
    Rcpp::stop("Output file name cannot be empty.");
  }
  filter_stream_in instream;
  open_fastx_in(instream, file);
  if (!instream) {
    Rcpp::stop("Error opening input file %s for reading", file);
  }

  bool gzipped = false;
  if (output.size() > 3 && output.substr(output.size() - 3) == ".gz") {
    gzipped = true;
  }

  filter_stream_out outstream;
  open_fastx_out(outstream, output, gzipped);
  if (!outstream) {
    Rcpp::stop("Error opening output file %s for writing", output);
  }

  // temporary storage for input records
  std::vector<std::string> pool;
  pool.reserve(denominator);

  std::string record;
  while (std::getline(instream, record) && !record.empty()){
    if (rename) {
      // Rename the read with a hexadecimal sequential number
      static std::size_t read_count = 0;
      std::stringstream ss;
      ss << "@" << std::hex << ++read_count;
      record = ss.str();
    }
    std::string line;
    std::getline(instream, line); // Sequence
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ record.");
    }
    record += "\n" + line;
    std::getline(instream, line); // Header2
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ record.");
    }
    record += "\n" + line;
    std::getline(instream, line); // Quality
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ record.");
    }
    record += "\n" + line;

    pool.push_back(record);

    if (pool.size() == denominator) {
      // Reservoir sampling: randomly select a record to keep
      Rcpp::IntegerVector idx = Rcpp::sample(denominator, denominator);
      for (int i = 0; i < numerator; ++i) {
        int selected = idx[i] - 1; // Convert to zero-based index
        outstream << pool[selected] << "\n";
      }
      pool.clear();
    }
  }
  if (pool.size() > 0) {
    // If there are remaining records, sample from them
    Rcpp::IntegerVector idx = Rcpp::sample(denominator, denominator);
    for (int i = 0; i < numerator; ++i) {
      int selected = idx[i] - 1; // Convert to zero-based index
      if (selected < pool.size()) {
        outstream << pool[selected] << "\n";
      }
    }
  }

  return output;
}

//' Take multiple samples from a FASTQ file
//'
//' This function takes multiple samples from a FASTQ file using the same
//' denomerator and seed, but different numerators.
//'
//' @param file (`character` string) the path to the input FASTQ file.
//' @param numerators (`integer` vector) a vector of numerators for the
//' fractions of reads to sample.
//' @param denominator (`integer` scalar) the denominator of the fraction of
//' reads to sample.
//' @param output (`character` vector) a vector of output file names. If an
//' element ends in ".gz", the output will be gzipped.
//'
//' @return a character vector of output file names
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector fastq_sample_multiple(
    const std::string& file,
    const Rcpp::IntegerVector& numerators,
    const int denominator,
    const Rcpp::CharacterVector& output,
    bool rename = false
) {
  if (numerators.size() != output.size()) {
    Rcpp::stop("The length of numerators and output must be the same.");
  }
  if (denominator <= 0) {
    Rcpp::stop("Denominator must be a positive integer.");
  }
  int max_numerator = 0;
  for (const auto& num : numerators) {
    if (num <= 0 || num > denominator) {
      Rcpp::stop("Each numerator must be a positive integer and less than or"
                   "equal to the denominator.");
    }
    if (num > max_numerator) {
      max_numerator = num;
    }
  }

  if (file.empty()) {
    Rcpp::stop("Input file name cannot be empty.");
  }

  filter_stream_in instream;
  open_fastx_in(instream, file);

  std::vector<filter_stream_out> out_streams(output.size());
  for (size_t i = 0; i < output.size(); ++i) {
    if (output[i].empty()) {
      Rcpp::stop("Output file name cannot be empty.");
    }
    std::string output_str(output[i]);
    open_fastx_out(out_streams[i], output_str);
    if (!out_streams[i]) {
      Rcpp::stop("Error opening output file %s for writing", output_str);
    }
  }

  // temporary storage for input records
  std::vector<std::string> pool;
  pool.reserve(denominator);

  std::string record;
  while (std::getline(instream, record) && !record.empty()) {
    if (rename) {
      // Rename the read with a hexadecimal sequential number
      static std::size_t read_count = 0;
      std::stringstream ss;
      ss << "@" << std::hex << ++read_count;
      record = ss.str();
    }
    std::string line;
    std::getline(instream, line); // Sequence
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ record.");
    }
    record += "\n" + line;
    std::getline(instream, line); // Header2
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ record.");
    }
    record += "\n" + line;
    std::getline(instream, line); // Quality
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ record.");
    }
    record += "\n" + line;

    pool.push_back(record);

    if (pool.size() == denominator) {
      // Reservoir sampling: randomly select records to keep
      Rcpp::IntegerVector idx = Rcpp::sample(denominator, denominator);
      for (int i = 0; i < max_numerator; ++i) {
        int selected = idx[i] - 1; // Convert to zero-based index
        for (int j = 0; j < numerators.size(); ++j) {
          if (i < numerators[j]) {
            out_streams[j] << pool[selected] << "\n";
          }
        }
      }
      pool.clear();
    }
  }
  if (pool.size() > 0) {
    // If there are remaining records, sample from them
    Rcpp::IntegerVector idx = Rcpp::sample(denominator, denominator);
    for (int i = 0; i < max_numerator; ++i) {
      int selected = idx[i] - 1; // Convert to zero-based index
      if (selected >= pool.size()) {
        continue; // Skip if selected index is out of bounds
      }
      for (int j = 0; j < numerators.size(); ++j) {
        if (i < numerators[j] && selected < pool.size()) {
          out_streams[j] << pool[selected] << "\n";
        }
      }
    }
  }

  return output;
}
