#include <Rcpp.h>

//' Remove unneeded fields from derep object(s)
//' @param derep (`list`) [`derep`][dada2::derep-class] object, or possibly
//' nested list of such objects
//' @return None. This function is called for its side effect, which is to
//' remove the `uniques` and `quals` fields from the input object
//' @export
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
