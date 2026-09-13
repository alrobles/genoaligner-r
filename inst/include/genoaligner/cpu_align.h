// genoaligner — CPU alignment core for the R package.
//
// This is the CPU backend that ships on CRAN. CRAN builds with a plain C++
// (g++) toolchain and never has a GPU, so this file is a self-contained C++17
// implementation of the two alignment algorithms the library exposes. It makes
// no reference to HIP or CUDA, so it compiles and runs everywhere.
//
// The GPU backend (ROCm/CUDA) of the sibling C++ library implements the SAME
// algorithms and agrees with this CPU core on the answers; wavefront
// acceleration and device speed are the GPU layer's job, not this file's.
// This core's job is CORRECTNESS on CPU.
//
// CIGAR CONVENTION (kept identical across both entry points so the package is
// self-consistent, and matching the C++ project's SW reference):
//     M  match       consumes a text base and a pattern base
//     X  mismatch    consumes a text base and a pattern base
//     I  insertion   consumes a TEXT (row/reference) base only
//     D  deletion    consumes a PATTERN (col/query) base only
// For full-length edit distance the orientation is a labeling of the same
// symmetric distance; for local SW it covers exactly the aligned span
// text[start_i..end_i] and pattern[start_j..end_j].
//
// WFA-derivation notice: the edit-distance DP computes the same Levenshtein
// distance the WFA formulation computes (sub=1, ins=1, del=1), validated
// against the trivial O(nm) rolling-row oracle below. The wavefront/WFA speed
// is not part of the CPU core; see the WFA attribution in the C++ library.
//
// The Smith-Waterman core is a port of the project's validated CPU reference
// (tests/sw/test_sw_trace.cpp, ref_sw_trace): full H/E/F matrices plus a
// value-comparison traceback walk, correct-beats-fast.

#ifndef GENOALIGNER_CPU_ALIGN_HPP
#define GENOALIGNER_CPU_ALIGN_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace genoaligner_cpu {

// SW op codes, matching the C++ project (include/genoaligner/backend/sw_kernel.hip).
constexpr int SW_OP_M = 0;
constexpr int SW_OP_X = 1;
constexpr int SW_OP_I = 2;   // consumes a text base
constexpr int SW_OP_D = 3;   // consumes a pattern base

// ---------------------------------------------------------------------------
// Oracle: classic Levenshtein, rolling row. Used ONLY to validate align_edit
// in tests (a second, obviously-correct implementation, like the C++ project's
// edit_distance_cpu). Sub=M=1, ins=1, del=1.
// ---------------------------------------------------------------------------
inline int edit_distance_oracle(const std::string& pattern, const std::string& text)
{
    const int m = (int)pattern.size(), n = (int)text.size();
    if (m == 0) return n;
    if (n == 0) return m;
    std::vector<int> prev(n + 1), cur(n + 1);
    for (int j = 0; j <= n; ++j) prev[j] = j;
    for (int i = 1; i <= m; ++i) {
        cur[0] = i;
        for (int j = 1; j <= n; ++j) {
            const int sub = prev[j - 1] + (pattern[(size_t)i - 1] == text[(size_t)j - 1] ? 0 : 1);
            cur[j] = std::min({sub, prev[j] + 1, cur[j - 1] + 1});
        }
        std::swap(prev, cur);
    }
    return prev[n];
}

// ---------------------------------------------------------------------------
// Edit distance (pattern vs text) with CIGAR, bounded by smax.
// ---------------------------------------------------------------------------
struct EditScore {
    int         score       = -1;   // edit distance; -1 if > smax (unresolved)
    bool        resolved    = false;
    std::string cigar;              // empty when unresolved or !with_cigar
    bool        rescore_ok  = false;
    bool        wellformed_ok = false;
};

inline EditScore align_edit(const std::string& pattern, const std::string& text,
                            int smax, bool with_cigar)
{
    EditScore r;
    const int m = (int)text.size(), n = (int)pattern.size();  // text = rows (m), pattern = cols (n)
    // smax out of [0,511] is a caller error -> unresolved (R layer rejects it).
    if (smax < 0 || smax > 511) return r;
    const int cap = smax + 1;

    // Degenerate empty inputs.
    if (m == 0 && n == 0) { r.score = 0; r.resolved = true; return r; }

    // Full DP with cap = smax+1 so we only distinguish <= smax from > smax,
    // and keep a direction byte per cell for the traceback.
    // text = rows (i in [0,m]), pattern = cols (j in [0,n]).
    std::vector<int>  dp((size_t)(m + 1) * (size_t)(n + 1));
    std::vector<uint8_t> dir((size_t)(m + 1) * (size_t)(n + 1), 0); // 0=diag(M/X),1=I(text),2=D(pattern)
    auto at = [&](std::vector<int>& M, int i, int j) -> int& { return M[(size_t)i * (size_t)(n + 1) + j]; };

    for (int j = 0; j <= n; ++j) at(dp, 0, j) = std::min(j, cap);       // no text rows: j pattern-moves (D)
    for (int i = 0; i <= m; ++i) at(dp, i, 0) = std::min(i, cap);       // no pattern cols: i text-moves (I)

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            const int diag = at(dp, i - 1, j - 1) + (text[(size_t)i - 1] == pattern[(size_t)j - 1] ? 0 : 1);
            const int up   = at(dp, i - 1, j) + 1;     // consume text base -> I
            const int left = at(dp, i, j - 1) + 1;     // consume pattern base -> D
            int best = diag, d = 0;
            if (up < best) { best = up; d = 1; }
            if (left < best) { best = left; d = 2; }
            at(dp, i, j) = std::min(best, cap);
            dir[(size_t)i * (size_t)(n + 1) + j] = (uint8_t)d;
        }
    }

    const int d = at(dp, m, n);
    if (d > smax) return r;                       // unresolved by design
    r.score = d;
    r.resolved = true;
    if (!with_cigar) return r;

    // Traceback.
    std::string rev;
    int i = m, j = n;
    while (i > 0 || j > 0) {
        const int d2 = (i > 0 && j > 0) ? (int)dir[(size_t)i * (size_t)(n + 1) + j] : 0;
        if (i > 0 && j > 0 && d2 == 0) {          // diagonal
            rev.push_back(text[(size_t)i - 1] == pattern[(size_t)j - 1] ? 'M' : 'X');
            --i; --j;
        } else if (i > 0 && (j == 0 || d2 == 1)) {  // consume text base -> I
            rev.push_back('I'); --i;
        } else {                                   // j > 0: consume pattern base -> D
            rev.push_back('D'); --j;
        }
    }
    std::reverse(rev.begin(), rev.end());
    r.cigar = rev;

    // Validation, computed in the library (same discipline as the C++ API).
    int ti = 0, pj = 0, rescore = 0; bool wf = true;
    for (char c : rev) {
        if (c == 'M' || c == 'X') {
            if (ti >= m || pj >= n) { wf = false; break; }
            const bool eq = text[(size_t)ti] == pattern[(size_t)pj];
            if ((c == 'M') != eq) wf = false;
            rescore += eq ? 0 : 1; ++ti; ++pj;
        } else if (c == 'I') { if (ti >= m) { wf = false; break; } ++ti; rescore += 1; }
        else { if (pj >= n) { wf = false; break; } ++pj; rescore += 1; }
    }
    wf = wf && ti == m && pj == n;
    r.wellformed_ok = wf;
    r.rescore_ok = (rescore == d);
    return r;
}

// ---------------------------------------------------------------------------
// Smith-Waterman (local, affine gaps) with CIGAR + coordinates.
// Port of the project's validated CPU reference (ref_sw_trace). text = rows (m),
// pattern = cols (n). scoring = {match, mismatch, gap_open, gap_extend};
// penalties are positive and subtracted.
// ---------------------------------------------------------------------------
struct SWParams { int match = 2, mismatch = -3, gap_open = 5, gap_extend = 2; };

struct SWResult {
    int         score       = -1;   // alignment score, >=0; -1 only if not computed
    bool        resolved    = false;
    std::string cigar;              // covers the aligned span; empty when score==0
    int         start_i = -1, start_j = -1, end_i = -1, end_j = -1;
    bool        rescore_ok = false, wellformed_ok = false;
    bool        too_large  = false; // matrix exceeded SW_MAX_TRACE_CELLS
};

// Largest direction table a pair may need (matches the C++ API limit).
constexpr long SW_MAX_TRACE_CELLS = 1L << 26;

inline SWResult align_sw(const std::string& text, const std::string& pattern, SWParams p,
                         bool with_cigar)
{
    SWResult r;
    const int m = (int)text.size(), n = (int)pattern.size();
    if (m == 0 || n == 0) { r.resolved = true; r.score = 0; return r; }
    if ((long)(m + 1) * (long)(n + 1) > SW_MAX_TRACE_CELLS) { r.too_large = true; r.resolved = false; return r; }

    const long W = n + 1L;
    std::vector<int> H((size_t)(m + 1) * W, 0), E(H.size(), 0), F(H.size(), 0);
    auto at = [&](std::vector<int>& M, int i, int j) -> int& { return M[(size_t)i * W + j]; };

    int best = 0, bi = 1, bj = 1;
    for (int i = 1; i <= m; ++i)
        for (int j = 1; j <= n; ++j) {
            const int s = (text[(size_t)i - 1] == pattern[(size_t)j - 1]) ? p.match : p.mismatch;
            int e = at(H, i, j - 1) - p.gap_open;
            { const int t = at(E, i, j - 1) - p.gap_extend; if (t > e) e = t; }
            int f = at(H, i - 1, j) - p.gap_open;
            { const int t = at(F, i - 1, j) - p.gap_extend; if (t > f) f = t; }
            int h = at(H, i - 1, j - 1) + s;
            if (e > h) h = e;
            if (f > h) h = f;
            if (h < 0) h = 0;
            at(H, i, j) = h; at(E, i, j) = e; at(F, i, j) = f;
            if (h > best) { best = h; bi = i; bj = j; }
        }
    r.score = best;
    r.resolved = true;
    if (best == 0) { r.start_i = r.start_j = r.end_i = r.end_j = -1; return r; }

    // Value-comparison traceback, ported verbatim from the validated reference.
    std::vector<int> codes;
    int i = bi, j = bj, state = 0;
    while (i > 0 && j > 0) {
        if (state == 0) {
            const int h = at(H, i, j);
            if (h == 0) break;
            const int s = (text[(size_t)i - 1] == pattern[(size_t)j - 1]) ? p.match : p.mismatch;
            if (h == at(H, i - 1, j - 1) + s) {
                codes.push_back(text[(size_t)i - 1] == pattern[(size_t)j - 1] ? SW_OP_M : SW_OP_X);
                --i; --j;
            } else if (h == at(E, i, j)) {
                codes.push_back(SW_OP_D);
                const bool ext = (at(E, i, j - 1) - p.gap_extend) > (at(H, i, j - 1) - p.gap_open);
                --j; state = ext ? 1 : 0;
            } else {  // h == F
                codes.push_back(SW_OP_I);
                const bool ext = (at(F, i - 1, j) - p.gap_extend) > (at(H, i - 1, j) - p.gap_open);
                --i; state = ext ? 2 : 0;
            }
        } else if (state == 1) {
            codes.push_back(SW_OP_D);
            const bool ext = (at(E, i, j - 1) - p.gap_extend) > (at(H, i, j - 1) - p.gap_open);
            --j; if (!ext) state = 0;
        } else {
            codes.push_back(SW_OP_I);
            const bool ext = (at(F, i - 1, j) - p.gap_extend) > (at(H, i - 1, j) - p.gap_open);
            --i; if (!ext) state = 0;
        }
    }
    r.start_i = i;   // text index where the aligned span starts
    r.start_j = j;   // pattern index
    r.end_i   = bi - 1;
    r.end_j   = bj - 1;
    if (!with_cigar) return r;

    // codes are reversed (element 0 = last op); emit the forward string.
    std::string cigar; cigar.resize(codes.size());
    for (int k = 0; k < (int)codes.size(); ++k) {
        const int c = codes[(size_t)codes.size() - 1 - (size_t)k];
        cigar[(size_t)k] = (c == SW_OP_M) ? 'M' : (c == SW_OP_X) ? 'X' : (c == SW_OP_I) ? 'I' : 'D';
    }
    r.cigar = cigar;

    // Well-formed + rescore over the aligned span. I consumes text, D consumes
    // pattern; a run of k same-direction gap ops = one gap of length k costing
    // gap_open + (k-1)*gap_extend.
    int ti = r.start_i, pj = r.start_j;
    int rescore = 0; bool wf = true;
    for (int k = 0; k < (int)cigar.size(); ++k) {
        const char c = cigar[(size_t)k];
        if (c == 'M' || c == 'X') {
            if (ti > r.end_i || pj > r.end_j) { wf = false; break; }
            const bool eq = text[(size_t)ti] == pattern[(size_t)pj];
            if ((c == 'M') != eq) wf = false;
            rescore += eq ? p.match : p.mismatch; ++ti; ++pj;
        } else {
            // gap run: count consecutive same-op runs
            char runc = c; int len = 0;
            while (k < (int)cigar.size() && cigar[(size_t)k] == runc) { ++len; ++k; }
            --k;
            rescore -= p.gap_open + (len - 1) * p.gap_extend;
            if (runc == 'I') { if (ti + len - 1 > r.end_i) { wf = false; break; } ti += len; }
            else             { if (pj + len - 1 > r.end_j) { wf = false; break; } pj += len; }
        }
    }
    wf = wf && ti == r.end_i + 1 && pj == r.end_j + 1;
    r.wellformed_ok = wf;
    r.rescore_ok = (rescore == r.score);
    return r;
}

} // namespace genoaligner_cpu

#endif // GENOALIGNER_CPU_ALIGN_HPP