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
Rcpp::CharacterVector fastq_sample_fraction(
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
Rcpp::CharacterVector fastq_sample_fraction_multiple(
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

//' Subsample a FASTQ file to a fixed number of reads
//'
//' This function subsamples a FASTQ file to a fixed number of reads.
//'
//' @param file (`character` string) the path to the input FASTQ file. May be gzipped.
//' @param number (`integer` scalar) the number of reads to sample.
//' @param output (`character` string) the path to the output file. If it ends in ".gz", the output will be gzipped.
//' @param rename (`logical` scalar) whether to rename the reads in the output;
//' if `TRUE`, the read names will be replaced with a hexadecimal sequential
//' number
//'
//' @return the output file name
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector fastq_sample_number(
    const std::string& file,
    const int number,
    const std::string& output,
    bool rename = false
) {

  int n_lines = 0;
  {
    filter_stream_in instream;
    open_fastx_in(instream, file);
    if (!instream) {
      Rcpp::stop("Error opening input file %s for reading", file);
    }
    std::string record;
    while (std::getline(instream, record) && !record.empty()) {
      n_lines++;
    }
  }
  if (n_lines % 4 != 0) {
    Rcpp::stop("Malformed FASTQ file: input file must have a multiple of 4 lines.");
  }

  int n_reads = n_lines / 4;
  filter_stream_in instream;
  open_fastx_in(instream, file);
  filter_stream_out outstream;
  open_fastx_out(outstream, output);
  Rcpp::IntegerVector sample_idx = Rcpp::sample(n_reads, n_reads);
  int i = 0;
  std::string header, sequence, header2, quality;
  while (std::getline(instream, header) && !header.empty()) {
    std::getline(instream, sequence);
    std::getline(instream, header2);
    std::getline(instream, quality);
    if (n_reads <= number || sample_idx[i++] <= number) {
      if (rename) {
        std::stringstream ss;
        ss << "@" << std::hex << i;
        header = ss.str();
      }
      outstream << header << "\n"
                << sequence << "\n"
                << header2 << "\n"
                << quality << "\n";
    }
  }

  return output;
}

//' Subsample a FASTQ file to multiple numbers of reads
//'
//' This function subsamples a FASTQ file to multiple numbers of reads.
//'
//' @param file (`character` string) the path to the input FASTQ file. May be gzipped.
//' @param numbers (`integer` vector) the numbers of reads to sample.
//' @param output (`character` vector) a vector of output file names. If an
//' element ends in ".gz", the output will be gzipped.
//' @param rename (`logical` scalar) whether to rename the reads in the output;
//' if `TRUE`, the read names will be replaced with a hexadecimal sequential
//' number.
//'
//' @return a character vector of output file names
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector fastq_sample_number_multiple(
    const std::string& file,
    const Rcpp::IntegerVector& numbers,
    const Rcpp::CharacterVector& output,
    bool rename = false
) {
  if (numbers.size() != output.size()) {
    Rcpp::stop("The length of numbers and output must be the same.");
  }
  if (file.empty()) {
    Rcpp::stop("Input file name cannot be empty.");
  }
  for (const auto& num : numbers) {
    if (num <= 0) {
      Rcpp::stop("Each number must be a positive integer.");
    }
  }

  int n_lines = 0;
  {
    filter_stream_in instream;
    open_fastx_in(instream, file);
    if (!instream) {
      Rcpp::stop("Error opening input file %s for reading", file);
    }
    std::string record;
    while (std::getline(instream, record) && !record.empty()) {
      n_lines++;
    }
  }
  if (n_lines % 4 != 0) {
    Rcpp::stop("Malformed FASTQ file: input file must have a multiple of 4 lines.");
  }

  int n_reads = n_lines / 4;
  Rcpp::IntegerVector sample_idx = Rcpp::sample(n_reads, n_reads);
  filter_stream_in instream;
  open_fastx_in(instream, file);
  std::vector<filter_stream_out> out_streams(output.size());
  for (int i = 0; i < output.size(); ++i) {
    if (output[i].empty()) {
      Rcpp::stop("Output file name cannot be empty.");
    }
    std::string output_str(output[i]);
    open_fastx_out(out_streams[i], output_str);
    if (!out_streams[i]) {
      Rcpp::stop("Error opening output file %s for writing", output_str);
    }
  }
  int i = 0;
  std::string header;
  while (std::getline(instream, header) && !header.empty()) {
    if (rename) {
      std::stringstream ss;
      ss << "@" << std::hex << (i + 1);
      header = ss.str();
    }
    std::string sequence, header2, quality;
    std::getline(instream, sequence);
    std::getline(instream, header2);
    std::getline(instream, quality);
    for (int j = 0; j < numbers.size(); ++j) {
      if (n_reads <= numbers[j] || sample_idx[i] <= numbers[j]) {
        out_streams[j] << header << "\n"
          << sequence << "\n"
          << header2 << "\n"
          << quality << "\n";
      }
    }
    i++;
  }
  return output;
}
