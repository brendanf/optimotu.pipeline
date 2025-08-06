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

//' @rdname fastq_sample_fraction
//' @param file_R1 (`character` string) the path to the R1 FASTQ file. May be gzipped.
//' @param file_R2 (`character` string) the path to the R2 FASTQ file. May be gzipped..
//' @param output_R1 (`character` string) the path to the R1 output file. If it ends in ".gz", the output will be gzipped.
//' @param output_R2 (`character` string) the path to the R2 output file. If it ends in ".gz", the output will be gzipped.
//' if `TRUE`, the read names will be replaced with a hexadecimal sequential
//' number.
//' @return a character vector of output file names
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector fastq_pair_sample_fraction(
    const std::string& file_R1,
    const std::string& file_R2,
    const int numerator,
    const int denominator,
    const std::string& output_R1,
    const std::string& output_R2,
    bool rename = false
) {


  if (numerator <= 0 || denominator <= 0 || numerator > denominator) {
    Rcpp::stop("Numerator and denominator must be positive integers, and"
                 " numerator must be less than or equal to denominator.");
  }
  if (file_R1.empty() || file_R2.empty()) {
    Rcpp::stop("Input file name cannot be empty.");
  }
  if (output_R1.empty() || output_R2.empty()) {
    Rcpp::stop("Output file name cannot be empty.");
  }
  filter_stream_in instream_R1, instream_R2;
  open_fastx_in(instream_R1, file_R1);
  if (!instream_R1) {
    Rcpp::stop("Error opening input file %s for reading", file_R1);
  }
  open_fastx_in(instream_R2, file_R2);
  if (!instream_R2) {
    Rcpp::stop("Error opening input file %s for reading", file_R2);
  }

  filter_stream_out outstream_R1, outstream_R2;
  open_fastx_out(outstream_R1, output_R1);
  if (!outstream_R1) {
    Rcpp::stop("Error opening output file %s for writing", output_R1);
  }
  open_fastx_out(outstream_R2, output_R2);
  if (!outstream_R2) {
    Rcpp::stop("Error opening output file %s for writing", output_R2);
  }

  // temporary storage for input records
  std::vector<std::pair<std::string, std::string>> pool;
  pool.reserve(denominator);

  std::string record_R1, record_R2;
  while (std::getline(instream_R1, record_R1) && !record_R1.empty()) {
    std::getline(instream_R2, record_R2);
    if (record_R2.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    if (rename) {
      // Rename the read with a hexadecimal sequential number
      static std::size_t read_count = 0;
      std::stringstream ss;
      ss << "@" << std::hex << ++read_count;
      record_R1 = ss.str();
      record_R2 = ss.str();
    }
    std::string line;
    std::getline(instream_R1, line); // Sequence
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    record_R1 += "\n" + line;
    std::getline(instream_R1, line); // Header2
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    record_R1 += "\n" + line;
    std::getline(instream_R1, line); // Quality
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    record_R1 += "\n" + line;

    std::getline(instream_R2, line); // Sequence
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    record_R2 += "\n" + line;
    std::getline(instream_R2, line); // Header2
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    record_R2 += "\n" + line;
    std::getline(instream_R2, line); // Quality
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    record_R2 += "\n" + line;
    pool.push_back(std::make_pair(record_R1, record_R2));

    if (pool.size() == denominator) {
      // Reservoir sampling: randomly select a record to keep
      Rcpp::IntegerVector idx = Rcpp::sample(denominator, denominator);
      for (int i = 0; i < numerator; ++i) {
        int selected = idx[i] - 1; // Convert to zero-based index
        outstream_R1 << pool[selected].first << "\n";
        outstream_R2 << pool[selected].second << "\n";
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
        outstream_R1 << pool[selected].first << "\n";
        outstream_R2 << pool[selected].second << "\n";
      }
    }
  }

  return Rcpp::CharacterVector::create(output_R1, output_R2);
}

//' @rdname fastq_sample_fraction_multiple
//' @param file_R1 (`character` string) the path to the R1 FASTQ file. May be gzipped.
//' @param file_R2 (`character` string) the path to the R2 FASTQ file. May be gzipped..
//' @param output_R1 (`character` vector) a vector of output file names. If an
//' element ends in ".gz", the output will be gzipped.
//' @param output_R2 (`character` vector) a vector of output file names. If an
//' element ends in ".gz", the output will be gzipped.
//' if `TRUE`, the read names will be replaced with a hexadecimal sequential
//' number.
//' @return a list of two character vectors of output file names
//' @export
// [[Rcpp::export]]
Rcpp::List fastq_pair_sample_fraction_multiple(
    const std::string& file_R1,
    const std::string& file_R2,
    const Rcpp::IntegerVector& numerators,
    const int denominator,
    const Rcpp::CharacterVector& output_R1,
    const Rcpp::CharacterVector& output_R2,
    bool rename = false
) {
  if (numerators.size() != output_R1.size() || numerators.size() != output_R2.size()) {
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

  if (file_R1.empty() || file_R2.empty()) {
    Rcpp::stop("Input file names cannot be empty.");
  }

  filter_stream_in instream_R1, instream_R2;
  open_fastx_in(instream_R1, file_R2);
  if (!instream_R1) {
    Rcpp::stop("Error opening input file %s for reading", file_R1);
  }
  open_fastx_in(instream_R2, file_R2);
  if (!instream_R2) {
    Rcpp::stop("Error opening input file %s for reading", file_R2);
  }

  std::vector<std::pair<filter_stream_out, filter_stream_out>> out_streams(output_R1.size());
  for (size_t i = 0; i < output_R1.size(); ++i) {
    if (output_R1[i].empty() || output_R2[i].empty()) {
      Rcpp::stop("Output file name cannot be empty.");
    }
    std::string output_str(output_R1[i]);
    open_fastx_out(out_streams[i].first, output_str);
    if (!out_streams[i].first) {
      Rcpp::stop("Error opening output file %s for writing", output_str);
    }
    output_str = output_R2[i];
    open_fastx_out(out_streams[i].second, output_str);
    if (!out_streams[i].second) {
      Rcpp::stop("Error opening output file %s for writing", output_str);
    }
  }

  // temporary storage for input records
  std::vector<std::pair<std::string, std::string>> pool;
  pool.reserve(denominator);

  std::string record_R1, record_R2;
  while (std::getline(instream_R1, record_R1) && !record_R1.empty()) {
    std::getline(instream_R2, record_R2);
    if (record_R2.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    if (rename) {
      // Rename the read with a hexadecimal sequential number
      static std::size_t read_count = 0;
      std::stringstream ss;
      ss << "@" << std::hex << ++read_count;
      record_R1 = ss.str();
      record_R2 = ss.str();
    }
    std::string line;
    std::getline(instream_R1, line); // Sequence
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    record_R1 += "\n" + line;
    std::getline(instream_R1, line); // Header2
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    record_R1 += "\n" + line;
    std::getline(instream_R1, line); // Quality
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    record_R1 += "\n" + line;

    std::getline(instream_R2, line); // Sequence
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    record_R2 += "\n" + line;
    std::getline(instream_R2, line); // Header2
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    record_R2 += "\n" + line;
    std::getline(instream_R2, line); // Quality
    if (line.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    record_R2 += "\n" + line;
    pool.push_back(std::make_pair(record_R1, record_R2));

    if (pool.size() == denominator) {
      // Reservoir sampling: randomly select records to keep
      Rcpp::IntegerVector idx = Rcpp::sample(denominator, denominator);
      for (int i = 0; i < max_numerator; ++i) {
        int selected = idx[i] - 1; // Convert to zero-based index
        for (int j = 0; j < numerators.size(); ++j) {
          if (i < numerators[j]) {
            out_streams[j].first << pool[selected].first << "\n";
            out_streams[j].second << pool[selected].second << "\n";
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
          out_streams[j].first << pool[selected].first << "\n";
          out_streams[j].second << pool[selected].second << "\n";
        }
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("R1") = output_R1, Rcpp::Named("R2") = output_R2);
}

//' @rdname fastq_sample_number
//' @param file_R1 (`character` string) the path to the R1 FASTQ file. May be gzipped.
//' @param file_R2 (`character` string) the path to the R2 FASTQ file. May be gzipped..
//' @param output_R1 (`character` string) the path to the R1 output file. If it ends in ".gz", the output will be gzipped.
//' @param output_R2 (`character` string) the path to the R2 output file. If it ends in ".gz", the output will be gzipped.
//' if `TRUE`, the read names will be replaced with a hexadecimal sequential
//' number.
//' @return a character vector of output file names
//' @export
// [[Rcpp::export]]
Rcpp::CharacterVector fastq_pair_sample_number(
    const std::string& file_R1,
    const std::string& file_R2,
    const int number,
    const std::string& output_R1,
    const std::string& output_R2,
    bool rename = false
) {

  int n_lines = 0;
  {
    filter_stream_in instream;
    open_fastx_in(instream, file_R1);
    if (!instream) {
      Rcpp::stop("Error opening input file %s for reading", file_R1);
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
  filter_stream_in instream_R1, instream_R2;
  open_fastx_in(instream_R1, file_R1);
  if (!instream_R1) {
    Rcpp::stop("Error opening input file %s for reading", file_R1);
  }
  open_fastx_in(instream_R2, file_R2);
  if (!instream_R2) {
    Rcpp::stop("Error opening input file %s for reading", file_R2);
  }
  filter_stream_out outstream_R1, outstream_R2;
  open_fastx_out(outstream_R1, output_R1);
  if (!outstream_R1) {
    Rcpp::stop("Error opening output file %s for writing", output_R1);
  }
  open_fastx_out(outstream_R2, output_R2);
  if (!outstream_R2) {
    Rcpp::stop("Error opening output file %s for writing", output_R2);
  }
  Rcpp::IntegerVector sample_idx = Rcpp::sample(n_reads, n_reads);
  int i = 0;
  std::string header_R1, header_R2, sequence_R1, sequence_R2, header2_R1, header2_R2, quality_R1, quality_R2;
  while (std::getline(instream_R1, header_R1) && !header_R1.empty()) {
    std::getline(instream_R2, header_R2);
    if (header_R2.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    std::getline(instream_R1, sequence_R1);
    if (sequence_R1.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    std::getline(instream_R1, header2_R1);
    if (header2_R1.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    std::getline(instream_R1, quality_R1);
    if (quality_R1.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    std::getline(instream_R2, sequence_R2);
    if (sequence_R2.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    std::getline(instream_R2, header2_R2);
    if (header2_R2.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    std::getline(instream_R2, quality_R2);
    if (quality_R2.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    if (n_reads <= number || sample_idx[i++] <= number) {
      if (rename) {
        std::stringstream ss;
        ss << "@" << std::hex << i;
        header_R1 = ss.str();
        header_R2 = ss.str();
      }
      outstream_R1 << header_R1 << "\n"
                << sequence_R1 << "\n"
                << header2_R1 << "\n"
                << quality_R1 << "\n";
      outstream_R2 << header_R2 << "\n"
                << sequence_R2 << "\n"
                << header2_R2 << "\n"
                << quality_R2 << "\n";
    }
  }

  return Rcpp::CharacterVector::create(output_R1, output_R2);
}

//' @rdname fastq_sample_number_multiple
//' @param file_R1 (`character` string) the path to the R1 FASTQ file. May be gzipped.
//' @param file_R2 (`character` string) the path to the R2 FASTQ file. May be gzipped..
//' @param output_R1 (`character` vector) a vector of output file names. If an
//' element ends in ".gz", the output will be gzipped.
//' @param output_R2 (`character` vector) a vector of output file names. If an
//' element ends in ".gz", the output will be gzipped.
//' @return a list of two character vectors of output file names
//' @export
// [[Rcpp::export]]
Rcpp::List fastq_pair_sample_number_multiple(
    const std::string& file_R1,
    const std::string& file_R2,
    const Rcpp::IntegerVector& numbers,
    const Rcpp::CharacterVector& output_R1,
    const Rcpp::CharacterVector& output_R2,
    bool rename = false
) {
  if (numbers.size() != output_R1.size() || numbers.size() != output_R2.size()) {
    Rcpp::stop("The length of numbers and output must be the same.");
  }
  if (file_R1.empty() || file_R2.empty()) {
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
    open_fastx_in(instream, file_R1);
    if (!instream) {
      Rcpp::stop("Error opening input file %s for reading", file_R1);
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
  filter_stream_in instream_R1, instream_R2;
  open_fastx_in(instream_R1, file_R1);
  if (!instream_R1) {
    Rcpp::stop("Error opening input file %s for reading", file_R1);
  }
  open_fastx_in(instream_R2, file_R2);
  if (!instream_R2) {
    Rcpp::stop("Error opening input file %s for reading", file_R2);
  }
  std::vector<std::pair<filter_stream_out, filter_stream_out>> out_streams(output_R1.size());
  for (int i = 0; i < output_R1.size(); ++i) {
    if (output_R1[i].empty() || output_R2[i].empty()) {
      Rcpp::stop("Output file name cannot be empty.");
    }
    std::string output_str(output_R1[i]);
    open_fastx_out(out_streams[i].first, output_str);
    if (!out_streams[i].first) {
      Rcpp::stop("Output file name cannot be empty.");
    }
    output_str = output_R2[i];
    open_fastx_out(out_streams[i].second, output_str);
    if (!out_streams[i].second) {
      Rcpp::stop("Error opening output file %s for writing", output_str);
    }
  }
  int i = 0;
  std::string header_R1, header_R2, sequence_R1, sequence_R2, header2_R1, header2_R2, quality_R1, quality_R2;
  while (std::getline(instream_R1, header_R1) && !header_R1.empty()) {
    std::getline(instream_R2, header_R2);
    if (header_R2.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    if (rename) {
      std::stringstream ss;
      ss << "@" << std::hex << (i + 1);
      header_R1 = ss.str();
      header_R2 = ss.str();
    }
    std::getline(instream_R1, sequence_R1);
    if (sequence_R1.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    std::getline(instream_R1, header2_R1);
    if (header2_R1.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    std::getline(instream_R1, quality_R1);
    if (quality_R1.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R1 record.");
    }
    std::getline(instream_R2, sequence_R2);
    if (sequence_R2.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    std::getline(instream_R2, header2_R2);
    if (header2_R2.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    std::getline(instream_R2, quality_R2);
    if (quality_R2.empty()) {
      Rcpp::stop("Unexpected end of file while reading FASTQ R2 record.");
    }
    for (int j = 0; j < numbers.size(); ++j) {
      if (n_reads <= numbers[j] || sample_idx[i] <= numbers[j]) {
        out_streams[j].first << header_R1 << "\n"
          << sequence_R1 << "\n"
          << header2_R1 << "\n"
          << quality_R1 << "\n";
        out_streams[j].second << header_R2 << "\n"
          << sequence_R2 << "\n"
          << header2_R2 << "\n"
          << quality_R2 << "\n";
      }
    }
    i++;
  }
  return Rcpp::List::create(Rcpp::Named("R1") = output_R1, Rcpp::Named("R2") = output_R2);
}

//' Sample a given number of reads from one or more FASTQ files
//'
//' This function samples the given number of reads from one or more FASTQ
//' files, based on a shuffle of the read indices. This guarantees reproducible
//' sampling when the same shuffle is reused. When multiple input files are
//' provided, the result will be as if the files were concatenated and then
//' sampled.
//' @param infile (`character`) the path(s) to the input FASTQ file(s). May be
//' gzipped.
//' @param outfile (`character` string) the path to the output FASTQ file. If it
//' ends in ".gz", the output will be gzipped.
//' @param n (`integer` scalar) the number of reads to sample.
//' @param sample (`integer` vector) a shuffle of the read indices.
//' @param rename (`logical`) if `TRUE`, the read names will be replaced with a hexadecimal sequential
//' number.
//' @return the output file name
//' @export
// [[Rcpp::export]]
std::string fastq_sample(
  const std::vector<std::string> & infile,
  const std::string & outfile,
  const int n,
  const Rcpp::IntegerVector sample,
  const bool rename = false
) {
  if (infile.empty() || outfile.empty()) {
    Rcpp::stop("Input and output file names cannot be empty.");
  }
  if (infile.size() == 1) {
    if (infile[0] == outfile) {
      if (sample.size() != n) {
        Rcpp::warning("Input (%s) and output (%s) file names are the same, but"
        " the length of sample (%d) is not equal to n (%d).",
        infile[0], outfile, sample.size(), n);
      }
      return outfile;
    }
  } else {
    for (const auto& file : infile) {
      if (file == outfile) {
        Rcpp::stop("Input and output file names cannot be the same when"
        " multiple input files are provided.");
      }
    }
  }
  if (n <= 0) {
    Rcpp::stop("n must be a positive integer.");
  }

  filter_stream_out outstream;
  open_fastx_out(outstream, outfile);
  if (!outstream) {
    Rcpp::stop("Error opening output file %s for writing", outfile);
  }

  filter_stream_in instream;
  int i = 0;
  for (const auto& file : infile) {
    filter_stream_in instream;
    open_fastx_in(instream, file);
    if (!instream) {
      Rcpp::stop("Error opening input file %s for reading", file);
    }
    while (i < sample.size()) {
      std::string record, line;
      if (sample[i] <= n) {
        std::getline(instream, record);
        if (record.empty()) {
          break;
        }
        if (rename) {
          std::stringstream ss;
          ss << "@" << std::hex << i;
          record = ss.str();
        }
        std::getline(instream, line); // sequence
        if (line.empty()) {
          Rcpp::stop("Unexpected end of file while reading FASTQ record.");
        }
        record += "\n" + line;
        std::getline(instream, line); // header2
        if (line.empty()) {
          Rcpp::stop("Unexpected end of file while reading FASTQ record.");
        }
        record += "\n" + line;
        std::getline(instream, line); // quality
        if (line.empty()) {
          Rcpp::stop("Unexpected end of file while reading FASTQ record.");
        }
        record += "\n" + line;
        outstream << record << "\n";
        if (!outstream) {
          Rcpp::stop("Error writing to output file %s", outfile);
        }
      } else {
        instream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if (instream.eof()) {
          break;
        }
        instream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        instream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        instream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        if (instream.eof()) {
          Rcpp::stop("Unexpected end of file while reading FASTQ record.");
        }
      }
      ++i;
    }
    if (i >= sample.size()) {
      break;
    }
  }
  return outfile;
}
