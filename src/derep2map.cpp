#include <Rcpp.h>

// [[Rcpp::export]]
void derep_map_only(Rcpp::List derep) {
  if (derep.inherits("derep")) {
    if (derep.containsElementNamed("uniques")) {
      derep["uniques"] = R_NilValue;
    }
    if (derep.containsElementNamed("quals")) {
      derep["quals"] = R_NilValue;
    }
  } else for (auto item : derep) {
    if (TYPEOF(item) == VECSXP) {
      derep_map_only(item);
    }
  }
}
