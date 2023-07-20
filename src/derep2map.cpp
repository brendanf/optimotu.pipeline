#include <Rcpp.h>

// [[Rcpp::export]]
void derep_map_only(Rcpp::List derep) {
  if (derep.inherits("derep")) {
    Rcpp::Rcout << "found a derep!";
    if (derep.containsElementNamed("uniques")) {
      Rcpp::Rcout << " removing uniques!";
      derep["uniques"] = R_NilValue;
    }
    if (derep.containsElementNamed("quals")) {
      Rcpp::Rcout << " removing quals!";
      derep["quals"] = R_NilValue;
    }
    Rcpp::Rcout << std::endl;
  } else for (auto item : derep) {
    if (TYPEOF(item) == VECSXP) {
      Rcpp::Rcout << "recursing...";
      derep_map_only(item);
    } else {
      Rcpp::Rcout << "element is not a list." << std::endl;
    }
  }
}
