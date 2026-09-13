// genoaligner — Rcpp glue. Thin wrappers: convert R scalars to the CPU core,
// pass results back as R lists. Batch/vectorisation is handled in R/; these
// handle ONE pair each so the R layer controls recycling and error reporting.
#include <Rcpp.h>
#include "genoaligner/cpu_align.h"

// [[Rcpp::export]]
Rcpp::List align_one(const std::string& pattern, const std::string& text,
                     int smax, bool with_cigar)
{
    Rcpp::List out;
    if (smax < 0 || smax > 511) {
        out["error"]   = "smax must be in [0, 511]";
        out["score"]   = NA_INTEGER;
        out["cigar"]   = NA_STRING;
        return out;
    }
    genoaligner_cpu::EditScore r =
        genoaligner_cpu::align_edit(pattern, text, smax, with_cigar);
    Rcpp::CharacterVector cg(1);
    if (r.resolved && with_cigar) { cg[0] = r.cigar; } else { cg[0] = NA_STRING; }
    out["score"]        = r.resolved ? r.score : NA_INTEGER;
    out["resolved"]     = r.resolved;
    out["cigar"]        = cg;
    out["rescore_ok"]   = r.rescore_ok;
    out["wellformed_ok"]= r.wellformed_ok;
    return out;
}

// [[Rcpp::export]]
Rcpp::List align_sw_one(const std::string& text, const std::string& pattern,
                        int match, int mismatch, int gap_open, int gap_extend,
                        bool with_cigar)
{
    genoaligner_cpu::SWParams p;
    p.match = match; p.mismatch = mismatch;
    p.gap_open = gap_open; p.gap_extend = gap_extend;
    genoaligner_cpu::SWResult r =
        genoaligner_cpu::align_sw(text, pattern, p, with_cigar);
    Rcpp::CharacterVector cg(1);
    if (r.resolved && !r.too_large && with_cigar) { cg[0] = r.cigar; }
    else { cg[0] = NA_STRING; }
    Rcpp::List out;
    out["score"]        = r.resolved && !r.too_large ? r.score : NA_INTEGER;
    out["resolved"]     = r.resolved;
    out["too_large"]    = r.too_large;
    out["start_i"]      = r.start_i;
    out["start_j"]      = r.start_j;
    out["end_i"]        = r.end_i;
    out["end_j"]        = r.end_j;
    out["cigar"]        = cg;
    out["rescore_ok"]   = r.rescore_ok;
    out["wellformed_ok"]= r.wellformed_ok;
    return out;
}