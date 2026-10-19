# Roadmap — codon-MSA improvement plan (post-2.0.0)

Grounded in `bench/out/sim_results{,auto,bbp}.tsv` and
`bench_results{,_fullscale}.tsv`. Ordered by expected impact.

## Where the measured gaps are

| regime | best ours | MAFFT | MACSE | decipher_seqs |
|---|---|---|---|---|
| d=0.1, pfs=0 | protein 0.994 | 0.999 | 0.983 | 0.997 |
| d=0.3, pfs=0 | protein 0.957 | 0.988 | 0.804 | 0.977 |
| d=0.1, pfs=0.2 | protein_auto 0.831 | 0.998 | 0.982 | 0.997 |
| d=0.3, pfs=0.2 | protein_auto 0.525 | 0.989 | 0.722 | 0.973 |
| fs recall (0.1/0.3) | protein_auto 0.31/0.36 | 0.46/0.48 | 0.74/0.31 | 0.44/0.59 |
| fs precision | ~0.01–0.02 | ~0.1–0.3 | ~0.37/0.06 | ~0.1–0.3 |

Two failure regimes: (1) pfs=0.2 collapse — architectural, frame-0
translation poisons everything downstream of an internal fs, and the
`len%%3` flag conflates terminal partials with internal events (on
real loci ~52% of rows are flagged → phase purity drops); (2) clean
d=0.3 SPS gap of ~0.03 vs MAFFT — shown *not* to be guide-tree,
`gap_open`, or `iter_refine` on the codon-token path.

## Phase 0 — make gaps measurable — DONE 2026-10-19 (diag_phase0.R)

**R1 DONE.** `simulate_codon.R` gained `--pcomp` (compensated fs pairs,
opposite sign 3–15 nt apart, net ≡0 mod3, shared `pair_id`), `--pterm`
(3' terminal trims = terminal partials, NOT events), `*.rowtruth.tsv`
(per-tip internal_fs/compensated/term_trim) and `*.tree.nwk`.

**R2 DONE. Answer: the guide tree is NOT the bottleneck.** `guide_tree=`
threads a phylo through `msa_align`/`msa_codon` (all backbones;
`msa_tree_cpp` shim + `phylo_to_merges`). True tree vs NJ over the
4-config grid: mean ΔSPS −0.003…+0.006, wins 13/24, loses 7 — noise.
The residual gap vs MAFFT is scoring/gap-model, not ordering.

**R3 DONE.** `msa_sp_score(x, mode)` exported (`msa_sp_score_cpp`;
codon mode tokenizes triplets via `cr_token3`, refuses non-mult-3
grids — subset `cofs` columns for refined MSAs). Model-SP vs true-SPS
rank agreement on NJ-vs-tree pairs: 0.917 — usable as a diagnostic,
not yet reliable enough alone for `select="sp"`.

**Flag audit (drives R4):** with only terminal trims (cfg C) the
`len%%3` flag false-fires on 11–17 rows/dataset — the real-data
confound reproduced in silico. With compensated pairs (cfg B) flag
recall for internal-fs rows is 0.11–0.64 — the invisible pairs are
quantified. Both failure modes live in `bench/out/diag_phase0.tsv`.

## Phase 1 — the frameshift regime (largest gap)

**R4. Guided re-tokenization flag. [M, shim on vendored machinery]**
`codon_encode_guided` is vendored; build leave-one-out guide bonus
from the pass-1 nt alignment (port the column-counting loop from
`genoaligner-devel/tools/genomsa.cpp`); `qc.partial` then counts
internal fs blocks + trailing block → `flag = fs_blocks >= 1 |
nchar%%3!=0`. Separates internal fs from terminal partials — the flag
is provably the binding constraint on compensated pairs (net
len%%3==0 invisible today) AND on ragged real data (~52% flagged now,
mostly terminal noise).
Fallback diagnostic: sliding-window stop-cluster score (pure R, weak
at d=0.3 but a useful QC cross-check).

**R5. Local-frame protein backbone — flagship. [M–L, package-side
composition]** Fix the root cause of the pfs=0.2 collapse:
`codon_encode_local`/`_guided` tile each row into in-frame codon
blocks and count fs blocks; translate only block codons; protein MSA;
`codon_backproject` unchanged; `codon_refine_sel_cpp` restores fs nts
as event columns (drift band `|s−3k|≤40` accommodates skipped nts).
Target: SPS ~0.9/0.7+ and fs_recall ~0.6+ at n=50, MACSE-class.
Risk: stop-avoiding encode DP can invent fs blocks at high divergence
(fs_enc=80 < 2×stop_pen=120) — mitigate via guided tiling / two-pass
(pass-1 frame-0 backbone → guide → pass-2 tiling) and `qc$fs_blocks`
provenance. Optional upstream `codon_encode_blocks()` returning fs
positions sharpens the flag.

**R6. fs-event precision package. [S–M, package-side]**
fs_prec is ~0.01–0.02 vs MACSE 0.2–0.5 — flat `codon_fs=30` is cheap
against a divergent profile. (a) per-row adaptive fs cost in
`sel_realign_row` (crp is built per row — scale by row divergence
evidence); (b) `rounds=2` arm; (c) two-threshold emission — rerun at
fs and 2×fs, emit only events surviving both.

## Phase 2 — clean-regime SPS

**R7. `iter_refine` retest on protein backbone. [S, one bench arm]**
Dead-end verdict was on codon-token; add `msa_codon_protein_iter2`
(iter_refine=2, fft_band=32) to `bench_sim.R`.

**R8. Engine-knob exposure + gap sweep. [S–M, shim]**
`gap_open=11/ge=1` are ClustalW protein presets, never swept here.
Shim builds `genomsa::Params` directly and calls `genomsa::msa_align`
(bypasses `MsaRequest`, which lacks knobs). Sweep
go∈{6,8,11,14}, ge∈{0.5,1,2}, gappy∈{0.9,0.95,0}.

**R9. Merge protein backbone + codon scoring. [M package; L upstream]**
(a) consensus-backbone shim: run both backbones, score each row's
placement under the codon-column profile, splice better per row, then
sel-refine. (b) upstream: AA-anchored codon DP — restrict the
codon-token progressive DP to bands around the protein backbone's
diagonals (`align_profiles_diag` exists). Only pursue (b) if (a) shows
the merge wins.

## Phase 3 — consolidation

**R10. Content gate in insertion-block merging. [M, package-side]**
`codon_refine_merge` co-columnizes same-length insertion blocks
blindly; at high divergence same-length blocks collide spuriously —
part of the fsP problem. Custom merge gated on nt agreement or guide
adjacency can live in the shim without touching vendored code.

**R11. Defaults + provenance + paper sync. [S]**
After R4–R6: consider `backbone="protein"` as default (it's also
faster on real loci: 0.24 vs 0.56 s median, 269 vs 430 s full-scale);
`qc$fs_blocks`/`qc$flag` columns; NEWS + paper update — the
"compensated-pair" limitation becomes a fixed item.

**R12. (Optional) codon-model simulation + axis C. [M]**
MG94-style subs + dN/dS FPR for reviewer credibility; defer unless
reviewers demand.

## Upstream vs package-side

Everything is package-side (R shims over already-vendored public
functions) except R9b (banded DP anchor plumbing → upstream PR +
re-vendor) and optional R5 `codon_encode_blocks()`. Per VENDORED.txt
discipline, `msa_rcpp.cpp` shims are package code.

## Documented dead ends (do not re-try)

- Global `refine=TRUE` on clean data (degrades SPS up to −0.26).
- Stop-count flagging (degenerates to global refine at d=0.3).
- Free/banded codon relocation in the phase DP.
- nt k-mer guide trees for codon-token path; `gap_open`/`iter_refine`
  on codon-token path; codon-matrix knob-twiddling (converges to
  protein mode anyway).

## Real-data finding driving R4

On the 31 real loci (`bench_real_variants.R`, 2026-10-19):
`msa_codon_protein` holds phase=1.000 everywhere (median 0.24 s), but
`+ refine="auto"` drops phase to median 0.524 — ~52% of real sequences
carry terminal partials, which the `nchar%%3` flag cannot distinguish
from internal frameshifts. The guided flag is the unblocking step.
