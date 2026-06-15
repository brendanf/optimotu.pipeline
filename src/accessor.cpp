#include "accessor.h"

Rcpp::IntegerVector integer_column(
    Rcpp::RObject df,
    const char* col_name,
    const char* df_name
) {
  if (!Rcpp::is<Rcpp::DataFrame>(df)) {
    Rcpp::stop("Object '%s' is not a data.frame.", df_name);
  }
  Rcpp::DataFrame df2 = Rcpp::as<Rcpp::DataFrame>(df);

  if (!df2.containsElementNamed(col_name)) {
    Rcpp::stop("Data frame '%s' has no column named '%s'.", df_name, col_name);
  }
  Rcpp::RObject col_raw = df2[col_name];

  if (!Rcpp::is<Rcpp::IntegerVector>(col_raw)) {
    Rcpp::stop("Column '%s' of data frame '%s' is not an integer vector.",
      col_name, df_name);
  }
  return Rcpp::as<Rcpp::IntegerVector>(col_raw);
}

Rcpp::NumericVector numeric_column(
    Rcpp::RObject df,
    const char* col_name,
    const char* df_name
) {
  if (!Rcpp::is<Rcpp::DataFrame>(df)) {
    Rcpp::stop("Object '%s' is not a data.frame.", df_name);
  }
  Rcpp::DataFrame df2 = Rcpp::as<Rcpp::DataFrame>(df);

  if (!df2.containsElementNamed(col_name)) {
    Rcpp::stop("Data frame '%s' has no column named '%s'.", df_name, col_name);
  }
  Rcpp::RObject col_raw = df2[col_name];

  if (!Rcpp::is<Rcpp::NumericVector>(col_raw)) {
    Rcpp::stop("Column '%s' of data frame '%s' is not a numeric vector.",
      col_name, df_name);
  }
  return Rcpp::as<Rcpp::NumericVector>(col_raw);
}

