// genoaligner — Rcpp glue for the public MSA API (host engine only).
//
// Thin wrapper, same discipline as align_rcpp.cpp: R scalars in, the C++
// library does the alignment AND the validation, an R list comes back.
// The engine entry is genoaligner::msa_align (src/api_msa.cpp, vendored
// from the C++ repo — see inst/VENDORED.txt), so the result contract the
// package exposes is exactly the one the C++ gates verify.
#include <Rcpp.h>
#include "genoaligner/api.hpp"

// [[Rcpp::export]]
Rcpp::List msa_run_cpp(const std::vector<std::string>& seqs,
                       const std::string& mode, int gc_def,
                       int codon_refine, bool local_frame)
{
    genoaligner::MsaRequest req;
    req.seqs             = seqs;
    req.gc_def           = gc_def;
    req.codon_refine     = codon_refine;
    req.codon_local_frame = local_frame;
    if      (mode == "dna")     req.mode = genoaligner::MsaMode::dna;
    else if (mode == "protein") req.mode = genoaligner::MsaMode::protein;
    else if (mode == "codon")   req.mode = genoaligner::MsaMode::codon;
    else {
        // R validates mode first; reaching this is a caller bug, refused
        // rather than guessed (same rule as the library's own checks).
        return Rcpp::List::create(
            Rcpp::_["status"]  = "invalid_argument",
            Rcpp::_["error"]   = "msa_run_cpp: unknown mode",
            Rcpp::_["aligned"] = Rcpp::CharacterVector(0),
            Rcpp::_["width"]   = 0,
            Rcpp::_["frame"]   = Rcpp::IntegerVector(0),
            Rcpp::_["stops"]   = Rcpp::IntegerVector(0),
            Rcpp::_["partial"] = Rcpp::IntegerVector(0));
    }

    genoaligner::MsaResult r = genoaligner::msa_align(req);

    const char* st = r.status == genoaligner::MsaResult::Status::ok        ? "ok"
                   : r.status == genoaligner::MsaResult::Status::invalid_argument
                                                                       ? "invalid_argument"
                   : r.status == genoaligner::MsaResult::Status::device_error
                                                                       ? "device_error"
                   :                                                     "error";
    Rcpp::String err = r.error ? Rcpp::String(r.error) : NA_STRING;

    const R_xlen_t n = (R_xlen_t)r.qc.size();
    Rcpp::IntegerVector frame(n), stops(n), partial(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        frame[i]   = r.qc[i].frame   < 0 ? NA_INTEGER : r.qc[i].frame;
        stops[i]   = r.qc[i].stops   < 0 ? NA_INTEGER : r.qc[i].stops;
        partial[i] = r.qc[i].partial < 0 ? NA_INTEGER : r.qc[i].partial;
    }

    return Rcpp::List::create(
        Rcpp::_["status"]  = st,
        Rcpp::_["error"]   = err,
        Rcpp::_["aligned"] = Rcpp::wrap(r.aligned),
        Rcpp::_["width"]   = r.width,
        Rcpp::_["frame"]   = frame,
        Rcpp::_["stops"]   = stops,
        Rcpp::_["partial"] = partial);
}
