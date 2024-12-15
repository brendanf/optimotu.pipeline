#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <array>
#include <cstdint>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string/predicate.hpp>

typedef boost::iostreams::filtering_stream<boost::iostreams::input>
filter_stream_in;
#define NO_LIMIT (std::numeric_limits<std::streamsize>::max())

// [[Rcpp::depends(BH)]]
// [[Rcpp::export]]
Rcpp::IntegerVector fastq_qual_bins(
    std::vector<std::string> fastq,
    int max_n = 1e9,
    int offset = 33
) {
   int i = 0;
   std::array<bool, 256> hits;
   hits.fill(false);

   for (const auto & f : fastq) {
      filter_stream_in instream;
      if (boost::ends_with(f, ".gz")) {
         instream.push(boost::iostreams::gzip_decompressor());
      }
      instream.push(
         boost::iostreams::file_source(f, std::ios_base::in | std::ios_base::binary)
      );
      while (!instream.eof()) {
         std::uint8_t c = 0;
         instream.ignore(NO_LIMIT, '\n');
         instream.ignore(NO_LIMIT, '\n');
         instream.ignore(NO_LIMIT, '\n');
         while ((c = instream.get()) != '\n') {
            if (instream.eof()) break;
            hits[c] = true;
         }
         ++i;
         if (i == max_n) break;
         Rcpp::checkUserInterrupt();
      }
      if (i == max_n) break;
   }
   Rcpp::IntegerVector out;
   for (int c = 0; c <=255; ++c) {
      if (hits[c]) out.push_back((int)c - offset);
   }
   return out;
}
