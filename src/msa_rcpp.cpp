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

// Sum-of-pairs objective of a finished MSA (genomsa::msa_sp_score), the
// same scorer iterative refinement maximizes internally. Exposed for
// diagnostics: e.g. does the model SP rank backbone variants the way
// true-SPS does, enabling "select = sp" model selection later.
// mode "codon" expects a pure codon grid (row width %3 == 0): each
// triplet is tokenized with cr_token3, "---" stays a gap, partial or
// ambiguous triplets become token 64 (scores 0 against everything) --
// the same semantics the stage-2 profile uses. For a refined MSA with
// event columns, subset the codon columns first (cofs).
// [[Rcpp::export]]
double msa_sp_score_cpp(const std::vector<std::string>& rows,
                        const std::string& mode, int gc_def)
{
    // The engine scorer trusts equal-width rows; enforce it here.
    const int W = rows.empty() ? 0 : (int)rows[0].size();
    for (const auto& r : rows)
        if ((int)r.size() != W)
            Rcpp::stop("msa_sp_score: ragged rows");
    if (mode == "protein")
        return genomsa::msa_sp_score(rows, genomsa::protein_params());
    if (mode == "codon") {
        if (W % 3 != 0)
            Rcpp::stop("msa_sp_score(codon): row width must be a "
                       "multiple of 3 -- subset codon columns (cofs) "
                       "of a refined MSA first");
        std::vector<std::string> tok(rows.size());
        for (size_t i = 0; i < rows.size(); ++i) {
            const char* p = rows[i].c_str();
            std::string& t = tok[i];
            t.reserve(W / 3);
            for (int j = 0; j < W; j += 3) {
                const char* q = p + j;
                if (q[0] == '-' && q[1] == '-' && q[2] == '-')
                    t.push_back('-');
                else
                    t.push_back((char)(128 + genoaligner::cr_token3(q)));
            }
        }
        return genomsa::msa_sp_score(tok, genomsa::codon_params(gc_def));
    }
    return genomsa::msa_sp_score(rows, genomsa::Params{});
}

// Guide-tree-override MSA (R2 diagnostic arm): the caller supplies the
// progressive guide as a post-order merge list over seq indices
// (mleft[k]/mright[k] = engine child ids of internal node n+k; children
// always carry smaller ids than their parent, so the last merge is the
// root join). Mode dispatch mirrors msa_run_cpp: dna/protein feed
// sequences directly; codon tokenizes (strict or local-frame), aligns
// tokens under the supplied tree, decodes, then runs codon_refine for
// the requested rounds (local_frame forces >= 1, same rule as the
// vendored run_codon). Params match params_for() exactly.
// [[Rcpp::export]]
Rcpp::List msa_tree_cpp(const std::vector<std::string>& seqs,
                        const std::string& mode, int gc_def,
                        int codon_refine, bool local_frame,
                        int iter_refine, int fft_band, int fft_lags,
                        double fft_min_rel,
                        const std::vector<int>& mleft,
                        const std::vector<int>& mright)
{
    const int n = (int)seqs.size();
    const int nint = (int)mleft.size();
    auto bad = [&](const char* msg) {
        return Rcpp::List::create(
            Rcpp::_["status"]  = "invalid_argument",
            Rcpp::_["error"]   = msg,
            Rcpp::_["aligned"] = Rcpp::CharacterVector(0),
            Rcpp::_["width"]   = 0,
            Rcpp::_["frame"]   = Rcpp::IntegerVector(0),
            Rcpp::_["stops"]   = Rcpp::IntegerVector(0),
            Rcpp::_["partial"] = Rcpp::IntegerVector(0));
    };
    if (nint != n - 1 || nint != (int)mright.size())
        return bad("msa_tree_cpp: a binary guide tree over n leaves "
                   "needs exactly n-1 merges");
    genomsa::Tree t;
    t.nodes.resize(n + nint);
    for (int k = 0; k < nint; ++k) {
        const int u = n + k;
        if (mleft[k] < 0 || mright[k] < 0 || mleft[k] >= u || mright[k] >= u)
            return bad("msa_tree_cpp: merge children must be existing "
                       "lower-id nodes (post-order)");
        t.nodes[u].left  = mleft[k];
        t.nodes[u].right = mright[k];
    }
    t.root = n + nint - 1;

    genomsa::Params P;
    std::vector<genomsa::CodonQc> cqc;
    const bool is_codon = mode == "codon";
    if      (mode == "protein") P = genomsa::protein_params();
    else if (is_codon) {
        P = genomsa::codon_params(gc_def);
        P.codon_refine      = codon_refine;
        P.codon_local_frame = local_frame ? 1 : 0;
    }
    else if (mode != "dna") return bad("msa_tree_cpp: unknown mode");
    P.iter_refine = iter_refine;
    P.fft_band    = fft_band;
    P.fft_lags    = fft_lags;
    P.fft_min_rel = (float)fft_min_rel;

    std::vector<std::string> msa;
    if (!is_codon) {
        msa = genomsa::msa_align_with_tree(seqs, t, P);
    } else {
        std::vector<std::string> tok =
            local_frame ? genomsa::codon_encode_local(seqs, P, &cqc)
                        : genomsa::codon_encode(seqs, gc_def, &cqc);
        msa = genomsa::codon_decode(
            genomsa::msa_align_with_tree(tok, t, P));
        int rounds = codon_refine;
        if (local_frame && rounds < 1) rounds = 1;
        if (rounds > 0)
            msa = genomsa::codon_refine(seqs, msa, P, rounds);
    }

    const R_xlen_t m = (R_xlen_t)n;
    Rcpp::IntegerVector frame(m, NA_INTEGER), stops(m, NA_INTEGER),
        partial(m, NA_INTEGER);
    for (R_xlen_t i = 0; i < m && i < (R_xlen_t)cqc.size(); ++i) {
        frame[i]   = cqc[i].frame;
        stops[i]   = cqc[i].stops;
        partial[i] = cqc[i].partial;
    }
    return Rcpp::List::create(
        Rcpp::_["status"]  = "ok",
        Rcpp::_["error"]   = NA_STRING,
        Rcpp::_["aligned"] = Rcpp::wrap(msa),
        Rcpp::_["width"]   = msa.empty() ? 0 : (int)msa[0].size(),
        Rcpp::_["frame"]   = frame,
        Rcpp::_["stops"]   = stops,
        Rcpp::_["partial"] = partial);
}

// ---- R4: guided re-tokenization flag ----------------------------------
// Port of the --guide-aln bonus construction from genoaligner-devel's
// tools/genomsa.cpp: per-column nt frequencies over non-gap rows; raw
// nt r of seq s maps to the column where its own row carries the r-th
// non-gap char, and bonus[r] = leave-one-out frequency - 0.5. Rows whose
// pass-1 representation under-covers the sequence (e.g. strict-encode
// lead drops) get an all-zero guide -- same "empty guide" semantics the
// vendored codon_encode_guided expects.
// [[Rcpp::export]]
Rcpp::List guide_bonus_cpp(const std::vector<std::string>& seqs,
                           const std::vector<std::string>& aln)
{
    const int W = aln.empty() ? 0 : (int)aln[0].size();
    const int n = (int)aln.size();
    std::vector<std::array<int,4>> cnt(W);
    std::vector<int> ngap(W, 0);
    for (const auto& row : aln)
        for (int c = 0; c < W && c < (int)row.size(); ++c) {
            int b = genoaligner::cr_nt4(row[c]);
            if (b >= 0) { ++cnt[c][b]; ++ngap[c]; }
        }
    Rcpp::List out(n);
    for (int s = 0; s < n; ++s) {
        const std::string& row = aln[s], &q = seqs[s];
        Rcpp::NumericVector b(q.size(), 0.0);
        int r = 0;
        for (int c = 0; c < W && c < (int)row.size() && r < (int)q.size();
             ++c) {
            if (row[c] == '-') continue;
            int b4 = genoaligner::cr_nt4(q[r]);
            if (b4 >= 0 && ngap[c] > 1) {
                float f = (float)(cnt[c][b4]
                        - (genoaligner::cr_nt4(row[c]) == b4 ? 1 : 0))
                        / (float)(ngap[c] - 1);
                b[r] = f - 0.5f;
            }
            ++r;
        }
        if (r == (int)q.size()) out[s] = b; else out[s] = Rcpp::NumericVector(0);
    }
    return out;
}

// Encode-DP replica emitting block provenance. The vendored
// codon_encode_guided already produces the token stream + qc (used here
// for parity checking); what it does NOT expose is WHERE the held-out
// frameshift/partial blocks sit, which is what the flag needs: internal
// fs blocks (codon blocks on both sides) are real frameshift evidence,
// while lead/trail blocks are the ragged-CDS terminal partials that the
// len%%3 flag wrongly counts. Tables mirror the vendored static ones
// (msa_ref.cpp) -- the parity test pins equivalence.
namespace {
const signed char SHIM_AA_STD[64] = {
    13,13,10,10, 15,15,15,15, 18,18,-1,-1,  4, 4,-1,17,
    10,10,10,10, 14,14,14,14,  8, 8, 5, 5,  1, 1, 1, 1,
     9, 9, 9,12, 16,16,16,16,  2, 2,11,11, 15,15, 1, 1,
    19,19,19,19,  0, 0, 0, 0,  3, 3, 6, 6,  7, 7, 7, 7};
const signed char SHIM_AA_MT[64] = {
    13,13,10,10, 15,15,15,15, 18,18,-1,-1,  4, 4,17,17,
    10,10,10,10, 14,14,14,14,  8, 8, 5, 5,  1, 1, 1, 1,
     9, 9,12,12, 16,16,16,16,  2, 2,11,11, 15,15,-1,-1,
    19,19,19,19,  0, 0, 0, 0,  3, 3, 6, 6,  7, 7, 7, 7};
}  // namespace

// [[Rcpp::export]]
Rcpp::List codon_encode_blocks_cpp(const std::vector<std::string>& seqs,
                                   int gc_def, double guide_w,
                                   const Rcpp::List& guide,
                                   double fs_enc = -1.0)
{
    const signed char* aa = gc_def == 2 ? SHIM_AA_MT : SHIM_AA_STD;
    const float NEG = -1e30f;
    genomsa::Params P = genomsa::codon_params(gc_def);
    const float fs_term = P.codon_fs_term,
                fs_cost = fs_enc >= 0 ? (float)fs_enc : P.codon_fs_enc,
                stop_pen = P.codon_stop_pen;
    const int n = (int)seqs.size();
    std::vector<std::vector<float>> gvec(n);
    for (int s = 0; s < n; ++s)
        if (s < guide.size())
            gvec[s] = Rcpp::as<std::vector<float>>(guide[s]);
    // Authoritative tokenization for parity assertions — only when the
    // DP parameters are the vendored defaults (fs_enc<0 keeps
    // P.codon_fs_enc); an override intentionally changes the tiling.
    std::vector<genomsa::CodonQc> qc;
    if (fs_enc < 0) {
        P.guide_w = (float)guide_w;
        (void)genomsa::codon_encode_guided(seqs, P, gvec, &qc);
    }

    Rcpp::IntegerVector n_fs(n), fs_nt(n), lead_nt(n), trail_nt(n),
        n_codons(n), vend_partial(n), parity(n);
    Rcpp::List fs_start(n), fs_len(n);
    for (int s = 0; s < n; ++s) {
        const std::string& q = seqs[s];
        const int L = (int)q.size();
        fs_start[s] = Rcpp::IntegerVector(0);
        fs_len[s]   = Rcpp::IntegerVector(0);
        if (L < 3) { trail_nt[s] = L; continue; }
        const std::vector<float>* g =
            (s < (int)gvec.size() && (int)gvec[s].size() >= L)
            ? &gvec[s] : nullptr;
        auto fsc = [&](int pos) {
            return (pos < 3 || pos > L - 3) ? fs_term : fs_cost; };
        std::vector<float> dp(L + 1, NEG);
        std::vector<uint8_t> tr(L + 1, 0);
        dp[0] = dp[1] = dp[2] = 0.f;
        for (int i = 3; i <= L; ++i) {
            float bv = NEG; uint8_t bt = 0;
            int a = genoaligner::cr_nt4(q[i-3]), b = genoaligner::cr_nt4(q[i-2]),
                c = genoaligner::cr_nt4(q[i-1]);
            float cs = (a < 0 || b < 0 || c < 0) ? 0.f
                       : (aa[a*16 + b*4 + c] < 0 ? -stop_pen : 1.f);
            if (g)
                cs += P.guide_w * ((*g)[i-3] + (*g)[i-2] + (*g)[i-1]);
            float pv = dp[i-3];
            if (pv > NEG/2 && pv + cs > bv) { bv = pv + cs; bt = 1; }
            for (int d = 1; d <= 2; ++d) {
                pv = dp[i - d];
                if (pv <= NEG/2) continue;
                float cand = pv - fsc(i - d);
                if (cand > bv) { bv = cand; bt = (uint8_t)(1 + d); }
            }
            dp[i] = bv; tr[i] = bt;
        }
        int ei = L; float eb = dp[L];
        for (int r = 1; r <= 2; ++r)
            if (dp[L-r] - fs_term > eb) { eb = dp[L-r] - fs_term; ei = L-r; }
        std::vector<int> blocks;
        int i = ei;
        while (i > 0) {
            uint8_t t = tr[i];
            if (t == 0) break;
            if (t == 1) { blocks.push_back(3); i -= 3; }
            else        { int d = t - 1; blocks.push_back(-d); i -= d; }
        }
        if (ei < L) blocks.push_back(-(L - ei));
        std::reverse(blocks.begin(), blocks.end());
        // Internal fs block: a held-out block with a codon block on
        // BOTH sides. Adjacent fs blocks can occur, so count codons,
        // not just neighboring indices. Held blocks sitting at an edge
        // (no codons on one side) are terminal raggedness: they extend
        // lead_nt / trail_nt even though they live inside `blocks`.
        const int nb = (int)blocks.size();
        std::vector<int> cod_before(nb), cod_after(nb);
        for (int k = 0, cb = 0; k < nb; ++k) {
            cod_before[k] = cb; if (blocks[k] == 3) ++cb;
        }
        for (int k = nb - 1, ca = 0; k >= 0; --k) {
            cod_after[k] = ca; if (blocks[k] == 3) ++ca;
        }
        std::vector<int> fsp, fsl;
        int pos = i, held = 0;
        lead_nt[s] = i;                       // leading drop (tr==0 break)
        for (int k = 0; k < nb; ++k) {
            int bl = blocks[k];
            if (bl == 3) { ++n_codons[s]; pos += 3; continue; }
            int d = -bl; ++held;
            if (cod_before[k] > 0 && cod_after[k] > 0) {
                ++n_fs[s]; fs_nt[s] += d;
                fsp.push_back(pos); fsl.push_back(d);
            } else if (cod_before[k] == 0) lead_nt[s] += d;
            else                         trail_nt[s] += d;
            pos += d;
        }
        fs_start[s] = Rcpp::wrap(fsp);
        fs_len[s]   = Rcpp::wrap(fsl);
        // vendored qc parity: its `partial` counts every held-out block
        // (internal fs + edge fs + trailing); the leading drop is not a
        // block there -- identical rule here.
        vend_partial[s] = s < (int)qc.size() ? qc[s].partial : -1;
        parity[s] = (vend_partial[s] < 0 || vend_partial[s] == held)
                    ? 1 : 0;
    }
    return Rcpp::List::create(
        Rcpp::_["n_fs"]         = n_fs,
        Rcpp::_["fs_nt"]        = fs_nt,
        Rcpp::_["lead_nt"]      = lead_nt,
        Rcpp::_["trail_nt"]     = trail_nt,
        Rcpp::_["n_codons"]     = n_codons,
        Rcpp::_["fs_start"]     = fs_start,
        Rcpp::_["fs_len"]       = fs_len,
        Rcpp::_["vend_partial"] = vend_partial,
        Rcpp::_["parity"]       = parity);
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
// offsets (cofs tracks real positions once event columns exist).
// Nucleotides living in insertion columns — between codon slots or at
// the edges — are re-emitted as ins blocks at the same boundary so an
// unflagged row never loses content.
genomsa::detail::Place sel_identity_place(const std::string& row,
                                        const std::vector<int>& cofs) {
    genomsa::detail::Place pl;
    const int C = (int)cofs.size();
    pl.fill.assign(C, "---");
    pl.ins.assign(C + 1, {});
    if (cofs[0] > 0)
        pl.ins[0].push_back(row.substr(0, cofs[0]));
    for (int j = 0; j < C; ++j) {
        pl.fill[j] = row.substr(cofs[j], 3);
        int next = j + 1 < C ? cofs[j + 1] : (int)row.size();
        int end = cofs[j] + 3;
        if (next > end)
            pl.ins[j + 1].push_back(row.substr(end, next - end));
    }
    return pl;
}
}  // namespace

// [[Rcpp::export]]
Rcpp::List codon_refine_sel_cpp(const std::vector<std::string>& seqs,
                                const std::vector<std::string>& aln,
                                const std::vector<bool>& flag,
                                int gc_def, int rounds, double fs = -1.0)
{
    using namespace genomsa;
    Params P = codon_params(gc_def);
    if (fs >= 0) P.codon_fs = (float)fs;
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
