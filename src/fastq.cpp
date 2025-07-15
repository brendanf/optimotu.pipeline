#include "fastq.h"

// disable warnings for use of auto_ptr and unnecessary parentheses in boost
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wparentheses"

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string/predicate.hpp>

#pragma GCC diagnostic pop

void open_fastx_in(filter_stream_in & instream, const std::string & infile) {
  if (boost::ends_with(infile, ".gz")) {
    instream.push(boost::iostreams::gzip_decompressor());
  }
  instream.push(
    boost::iostreams::file_source(infile, std::ios_base::in | std::ios_base::binary)
  );
  if (!instream) {
    Rcpp::stop("Failed to open input file " + infile);
  }
}

void open_fastx_out(filter_stream_out & outstream, const std::string & outfile, bool compress) {
  if (compress) {
    outstream.push(boost::iostreams::gzip_compressor());
  }
  outstream.push(
    boost::iostreams::file_sink(outfile, std::ios_base::out | std::ios_base::binary)
  );
  if (!outstream) {
    Rcpp::stop("Failed to open output file " + outfile);
  }
}

void open_fastx_out(filter_stream_out & outstream, const std::string & outfile) {
  if (boost::ends_with(outfile, ".gz")) {
    open_fastx_out(outstream, outfile, true);
  } else {
    open_fastx_out(outstream, outfile, false);
  }
}
