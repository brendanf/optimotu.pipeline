#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <deque>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string/predicate.hpp>

typedef boost::iostreams::filtering_stream<boost::iostreams::input>
filter_stream_in;
typedef std::string string;

#define NO_LIMIT (std::numeric_limits<std::streamsize>::max())

// [[Rcpp::depends(BH)]]
// [[Rcpp::export]]
Rcpp::CharacterVector fastq_names(std::string x) {
  std::deque<string> names;
  filter_stream_in instream;
  if (boost::ends_with(x, ".gz")) {
    instream.push(boost::iostreams::gzip_decompressor());
  }
  instream.push(
    boost::iostreams::file_source(x, std::ios_base::in | std::ios_base::binary)
  );
  while (instream.get()) {
    string line;
    instream >> line;
    if (line.size() == 0) break;
    names.push_back(line);
    instream.ignore(NO_LIMIT, '\n');
    instream.ignore(NO_LIMIT, '\n');
    instream.ignore(NO_LIMIT, '\n');
    instream.ignore(NO_LIMIT, '\n');
  }
  return Rcpp::wrap(names);
}

bool read_fastq_record(std::istream &stream, std::string &s) {
  if (stream.get() != '@') return false;
  stream >> s;
  stream.ignore(NO_LIMIT, '\n');
  stream.ignore(NO_LIMIT, '\n');
  stream.ignore(NO_LIMIT, '\n');
  stream.ignore(NO_LIMIT, '\n');
  return (bool)stream;
}

// [[Rcpp::export]]
Rcpp::DataFrame fastq_stage_map(std::string raw, std::vector<std::string> stages) {
  if (stages.size() > 8)
    Rcpp::stop("only up to 8 processing stages are allowed in fastq_map");

  int i = 0;
  filter_stream_in raw_in;
  std::vector<std::unique_ptr<filter_stream_in>> stage_in;
  if (boost::ends_with(raw, ".gz")) {
    raw_in.push(boost::iostreams::gzip_decompressor());
  }
  raw_in.push(boost::iostreams::file_source(raw, std::ios_base::binary));

  for (auto s : stages) {
    auto s_in = std::make_unique<filter_stream_in>();
    if (boost::ends_with(s, ".gz")) {
      s_in->push(boost::iostreams::gzip_decompressor());
    }
    s_in->push(boost::iostreams::file_source(s, std::ios_base::binary));
    stage_in.push_back(std::move(s_in));
  }

  std::string raw_i;
  std::vector<std::string> stage_i(stages.size());
  for (std::size_t si = 0; si < stages.size(); ++si) {
    read_fastq_record(*stage_in[si], stage_i[si]);
  }

  std::deque<unsigned char> outflags;
  std::deque<int> name_index;
  std::deque<int> count_index;
  bool use_name = true;

  while (read_fastq_record(raw_in, raw_i)) {
    char of = 0;
    for (std::size_t si = 0; si < stages.size(); ++si) {
      if (stage_in[si] && stage_i[si] == raw_i) {
        of |= (unsigned char)0x01 << si;
        read_fastq_record(*stage_in[si], stage_i[si]);
      }
    }
    outflags.push_back(of);

    count_index.push_back(++i);

    if (use_name) {
      std::stringstream ss(raw_i);
      int j;
      if (ss >> std::hex >> j) {
        name_index.push_back(j);
      } else {
        name_index.clear();
        use_name = false;
      }
    }
  }

  Rcpp::RawVector flags = Rcpp::wrap(outflags);
  Rcpp::IntegerVector index = Rcpp::wrap(use_name ? name_index : count_index);
  return Rcpp::DataFrame::create(
    Rcpp::Named("read_idx") = index,
    Rcpp::Named("flags") = flags
  );
}

// [[Rcpp::export]]
Rcpp::DataFrame fastq_stage_map2(std::string raw, std::string stage1, std::string stage2) {
  std::vector<std::string> stages{stage1, stage2};
  return fastq_stage_map(raw, stages);
}
