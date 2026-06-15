#ifndef OPTIMOTU_PIPELINE_ACCESSOR_H
#define OPTIMOTU_PIPELINE_ACCESSOR_H
#include <Rcpp.h>
#include <string>

Rcpp::IntegerVector integer_column(
    Rcpp::RObject df,
    const char* col_name,
    const char* df_name
);

Rcpp::NumericVector numeric_column(
    Rcpp::RObject df,
    const char* col_name,
    const char* df_name
);

#endif
