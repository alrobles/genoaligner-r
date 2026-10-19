// genoaligner — Rcpp glue for the public MSA API (host engine only).
//
// Thin wrapper, same discipline as align_rcpp.cpp: R scalars in, the C++
// library does the alignment AND the validation, an R list comes back.
// The engine entry is genoaligner::msa_align (src/api_msa.cpp, vendored
// from the C++ repo — see inst/VENDORED.txt), so the result contract the
// package exposes is exactly the one the C++ gates verify.
#include <Rcpp.h>
#include "genoaligner/api.hpp"
#include "genoaligner/msa/msa.hpp"
#include <genoaligner/backend/codon_refine_kernel_impl.hip>

// [[Rcpp::export]]
Rcpp::List msa_run_cpp(const std::vector<std::string>& seqs,
                       const std::string& mode, int gc_def,
                       int codon_refine, bool local_frame,
                       int iter_refine, int fft_band, int fft_lags,
                       double fft_min_rel)
{
    genoaligner::MsaRequest req;
    req.seqs             = seqs;
    req.gc_def           = gc_def;
    req.codon_refine     = codon_refine;
    req.codon_local_frame = local_frame;
    req.iter_refine      = iter_refine;
    req.fft_band         = fft_band;
    req.fft_lags         = fft_lags;
    req.fft_min_rel      = (float)fft_min_rel;
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

// bench-only shim: run the stage-2 codon phase DP over an arbitrary
// codon-grid alignment (e.g. a protein back-projection skeleton).
// [[Rcpp::export]]
Rcpp::CharacterVector codon_refine_cpp(const std::vector<std::string>& seqs,
                                       const std::vector<std::string>& aln,
                                       int gc_def, int rounds)
{
    genomsa::Params P = genomsa::codon_params(gc_def);
    return Rcpp::wrap(genomsa::codon_refine(seqs, aln, P, rounds));
}

// Selective stage-2 refinement: the per-row phase DP only runs for rows
// flagged true; every other row contributes its input placement verbatim
// (its fills in the codon columns, no event columns). The profile still
// sees the full alignment, so flagged rows realign against the same
// context as a global codon_refine -- only event ownership is gated.
// Body mirrors genomsa::codon_refine over public header pieces
// (detail::codon_prof_build / codon_phase_dp / codon_trace_place /
// codon_refine_merge); the vendored engine is untouched.
namespace {
genomsa::detail::Place sel_realign_row(
        const std::string& seq, const genomsa::detail::CodonProf& cp,
        int self, const genomsa::Params& P,
        const std::vector<float>& dot) {
    using namespace genomsa::detail;
    const int L = (int)seq.size(), C = cp.C;
    const std::vector<int>& stok = cp.tok[self];
    std::vector<int> pi;
    std::vector<int8_t> stok8(C);
    for (int j = 0; j < C; ++j) {
        stok8[j] = (int8_t)stok[j];
        if (stok[j] >= 0) pi.push_back(j);
    }
    const int K = (int)pi.size();
    if (L == 0 || C == 0 || cp.nseq <= 1 || K == 0)
        return codon_trace_place(seq, cp, self, P, pi, nullptr, 0, -2);
    genoaligner::CodonRefineParams crp{
        P.match, P.ts, P.tv, P.gap_open, P.psgp ? 1 : 0,
        P.psgp_scale, P.psgp_min_open,
        P.codon_fs, P.codon_fs_term,
        (int)P.codon_refine_band, cp.nseq, C};
    const int ND = 2 * crp.band + 1;
    std::vector<uint8_t> tr((size_t)K * ND, 0);
    std::vector<float>   dpb(2 * ND);
    const float* ncnt = &cp.ncnt[0][0][0];
    int bes = -1;
    genoaligner::codon_phase_dp(seq.c_str(), L, stok8.data(), pi.data(), K,
                                dot.data(), cp.ocnt.data(), ncnt,
                                P.sub.data(), crp, tr.data(), dpb.data(),
                                &bes);
    return codon_trace_place(seq, cp, self, P, pi, tr.data(), ND, bes);
}

// Identity placement: the row's input fills at the current codon
// offsets (cofs tracks real positions once event columns exist),
// no insertions.
genomsa::detail::Place sel_identity_place(const std::string& row,
                                        const std::vector<int>& cofs) {
    genomsa::detail::Place pl;
    const int C = (int)cofs.size();
    pl.fill.assign(C, "---");
    pl.ins.assign(C + 1, {});
    for (int j = 0; j < C; ++j) pl.fill[j] = row.substr(cofs[j], 3);
    return pl;
}
}  // namespace

// [[Rcpp::export]]
Rcpp::List codon_refine_sel_cpp(const std::vector<std::string>& seqs,
                                const std::vector<std::string>& aln,
                                const std::vector<bool>& flag,
                                int gc_def, int rounds)
{
    using namespace genomsa;
    Params P = codon_params(gc_def);
    std::vector<std::string> cur = aln;
    int C = (int)cur[0].size() / 3;
    std::vector<int> cofs(C);
    for (int j = 0; j < C; ++j) cofs[j] = 3 * j;
    for (int r = 0; r < rounds; ++r) {
        detail::CodonProf cp = detail::codon_prof_build(cur, cofs);
        std::vector<float> dot((size_t)65 * cp.C, 0.f);
        for (int t = 0; t < 65; ++t)
            for (int j = 0; j < cp.C; ++j) {
                float s = 0;
                for (int b = 0; b < 65; ++b)
                    if (cp.ccnt[j][b]) s += cp.ccnt[j][b] * P.sub[t * 65 + b];
                dot[(size_t)t * cp.C + j] = s;
            }
        std::vector<detail::Place> pls(cp.nseq);
        for (int s = 0; s < cp.nseq; ++s)
            pls[s] = flag[s] ? sel_realign_row(seqs[s], cp, s, P, dot)
                             : sel_identity_place(cur[s], cofs);
        cur = detail::codon_refine_merge(pls, cp.C, &cofs);
        C = (int)cofs.size();
    }
    return Rcpp::List::create(Rcpp::_["aln"]  = Rcpp::wrap(cur),
                              Rcpp::_["cofs"] = Rcpp::wrap(cofs));
}
