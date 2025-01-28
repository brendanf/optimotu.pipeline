#include "fastq.h"

// Internal function to read a FASTQ record
// @param stream (`&istream`) input stream
// @param s (`&string`) location to write header from next FASTQ record
// @return (`bool`) whether a record was read
bool read_fastq_record(std::istream &stream, std::string &s) {
  if (stream.get() != '@') return false;
  stream >> s;
  stream.ignore(NO_LIMIT, '\n');
  stream.ignore(NO_LIMIT, '\n');
  stream.ignore(NO_LIMIT, '\n');
  stream.ignore(NO_LIMIT, '\n');
  return (bool)stream;
}

//' Flag the fate of reads from a FASTQ file through processing stages
//' @param raw (`character`) FASTQ file path, optionally gzipped. The raw reads.
//' @param stages (`character`) FASTQ file paths, optionally gzipped. The
//' processed reads after multiple ordered stages.
//' @return a `data.frame` with columns `read_idx` (`integer`) and `flags`
//' (`raw`) indicating which processing stages each read was found in.
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame fastq_stage_flags(std::string raw, std::vector<std::string> stages) {
  if (stages.size() > 8)
    Rcpp::stop("only up to 8 processing stages are allowed in fastq_stage_flags");

  int i = 0;
  filter_stream_in raw_in;
  open_fastx_in(raw_in, raw);

  std::vector<filter_stream_in> stage_in(stages.size());
  for (size_t i = 0; i < stages.size(); ++i) {
    open_fastx_in(stage_in[i], stages[i]);
  }

  std::string raw_i;
  std::vector<std::string> stage_i(stages.size());
  for (std::size_t si = 0; si < stages.size(); ++si) {
    read_fastq_record(stage_in[si], stage_i[si]);
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
        read_fastq_record(stage_in[si], stage_i[si]);
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

//' @rdname fastq_stage_flags
//' @param stage1 (`character`) FASTQ file path, optionally gzipped, the first
//' processed stage
//' @param stage2 (`character`) FASTQ file path, optionally gzipped, the second
//' processed stage
//' @export
 // [[Rcpp::export]]
 Rcpp::DataFrame fastq_stage_flag2(std::string raw, std::string stage1, std::string stage2) {
  std::vector<std::string> stages{stage1, stage2};
  return fastq_stage_flags(raw, stages);
}

//' Map reads from a FASTQ file to processing stages
//' @param raw (`character`) FASTQ file path, optionally gzipped. The raw reads.
//' @param stages (`character`) FASTQ file paths, optionally gzipped. The
//' processed reads after multiple ordered stages.
//' @return a `data.frame` with columns `read_idx` (`integer`) and `flags`
//' (`raw`) indicating which processing stages each read was found in.
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame fastq_stage_map(std::string raw, Rcpp::CharacterVector stages) {
  std::vector<std::string> stages_v = Rcpp::as<std::vector<std::string>>(stages);

  int i = 0;
  filter_stream_in raw_stream;
  open_fastx_in(raw_stream, raw);

  std::vector<filter_stream_in> stage_stream(stages.size());
  for (int i = 0; i < stages.size(); ++i) {
    open_fastx_in(stage_stream[i], stages_v[i]);
  }

  std::string raw_name;
  std::vector<std::string> stage_name(stages_v.size());
  std::vector<int> stage_i(stages.size(), 1);
  for (int si = 0; si < stages.size(); ++si) {
    read_fastq_record(stage_stream[si], stage_name[si]);
  }

  std::deque<unsigned char> outflags;
  std::deque<int> name_index;
  std::deque<int> count_index;
  std::vector<std::deque<int>> stage_index(stages.size());
  bool use_name = true;

  while (read_fastq_record(raw_stream, raw_name)) {
    for (int si = 0; si < stages.size(); ++si) {
      if (stage_stream[si] && stage_name[si] == raw_name) {
        stage_index[si].push_back(stage_i[si]++);
        read_fastq_record(stage_stream[si], stage_name[si]);
      } else {
        stage_index[si].push_back(NA_INTEGER);
      }
    }

    count_index.push_back(++i);

    if (use_name) {
      std::stringstream ss(raw_name);
      int j;
      if (ss >> std::hex >> j) {
        name_index.push_back(j);
      } else {
        name_index.clear();
        use_name = false;
      }
    }
  }

  Rcpp::IntegerVector index = Rcpp::wrap(use_name ? name_index : count_index);
  Rcpp::DataFrame output = Rcpp::DataFrame::create(
    Rcpp::Named("raw_idx") = index
  );
  if (stages.hasAttribute("names")) {
    std::vector<std::string> stage_names = stages.names();
    for (int si = 0; si < stages.size(); ++si) {
      output.push_back(Rcpp::wrap(stage_index[si]), stage_names[si] + "_idx");
    }
  } else {
    for (int si = 0; si < stages.size(); ++si) {
      std::string name = "stage" + std::to_string(si + 1) + "_idx";
      output.push_back(Rcpp::wrap(stage_index[si]), name);
    }
  }
  return output;
}

//' Map the fate of individual reads through trimming and filtering steps
//' @param fq_raw (`character`) raw FASTQ file path, optionally gzipped
//' @param fq_trim (`character`) trimmed FASTQ file path, optionally gzipped
//' @param fq_filt (`character`) filtered FASTQ file path, optionally gzipped
//' @return `data.frame` with columns:
//' - `raw_idx` (integer) index of the sequence in fastq_raw; or if sequence names
//' in the fastq files are hex encoded integers (e.g., during subsampling to
//' indicate the original index) then `seq_id` comes from the sequence names.
//' - `trim_idx` (integer) index of the sequence in fastq_trim, or NA if absent
//' - `filt_idx` (integer) index of the sequence in fastq_filt, or NA if absent
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame fastq_seq_map(std::string fq_raw, std::string fq_trim, std::string fq_filt) {
  Rcpp::CharacterVector stages;
  stages["trim"] = fq_trim;
  stages["filt"] = fq_filt;
  return fastq_stage_map(fq_raw, stages);
}
