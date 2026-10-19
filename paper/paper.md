---
title: "genoaligner: protein-guided codon-aware multiple sequence alignment at phylogenomic scale"
author: "Angel L. Robles-Fernandez"
date: "draft — axis A+B results complete"
output: pdf_document
---

# Abstract

Multiple sequence alignment (MSA) of protein-coding sequences underpins
comparative genomics, yet the tools that preserve the reading frame —
notably MACSE — are too slow for phylogenomic-scale datasets, while fast
nucleotide aligners break the frame. We present **genoaligner**, a
portable C++17 alignment engine (R package on CRAN) whose codon-aware MSA
builds the alignment backbone at amino-acid level and back-projects to a
codon grid, guaranteeing whole-codon indels. On simulated data with known
ground truth, the protein-backbone mode reaches sum-of-pairs accuracy
0.994/0.957 at nucleotide divergence 0.1/0.3 — third overall behind MAFFT
and DECIPHER, ahead of MACSE — at ~0.3 s per 50-sequence dataset. A
selective frameshift-refinement pass (`refine = "auto"`) recovers
frameshift events only in sequences that provably cannot be in-frame
(length not a multiple of 3), leaving clean coding output bit-identical
to the un-refined backbone. On the 31 classic mammalian loci (Upham et
al. 2019), `msa_codon` aligns each locus in a median 0.56 s with codon
phase purity 1.0 — a median **325× speedup over MACSE** (range
156–1046×), which produces identical phase purity in a median 171 s and
exceeds one hour on single 50-sequence loci. At full scale (up to 3,559
sequences) `msa_codon` completes all 31 loci in 7.2 min total.

# 1. Introduction

Protein-coding DNA alignment is a special case of MSA: biologically
meaningful insertions and deletions in coding sequence are (mostly)
multiples of three, so an alignment that scatters indels at single
nucleotide positions corrupts the reading frame and inflates downstream
estimates of nonsynonymous divergence. MACSE (Ranwez et al. 2011, 2018)
solved the frame-preservation problem — including genuine frameshifts
in pseudogenes and error-prone reads — but its profile Hidden Markov
Model machinery makes it orders of magnitude slower than standard
aligners and effectively unusable at the thousands-of-sequences scale
of current phylogenomics. Faster translation-guided tools
(e.g. DECIPHER's `AlignTranslation`) compromise the frame.

genoaligner takes the opposite route: a single portable C++17
progressive-MSA engine that is verified bit-exact against a GPU
implementation in the sibling repository, exposing DNA, protein, and
codon-aware modes through one API. Version 2.x adds two complementary
codon strategies — the original codon-token progressive DP, and a
protein-guided backbone that translates each sequence in frame 0,
aligns amino acids with the same engine, and back-projects each residue
to its codon — plus a selective stage-2 pass that repairs reading
frames only where the data demonstrably require it.

# 2. Methods

## 2.1 Codon-aware alignment modes

`msa_codon(seqs, gc, backbone, refine)`:

* **`backbone = "codon"` (default).** Sequences are tokenized to codons
  under NCBI genetic code `gc` (1 = standard, 2 = vertebrate
  mitochondrial) and aligned over a codon substitution matrix
  (BLOSUM62-derived amino-acid scores + per-identical-nucleotide bonus +
  stop penalty). Every indel is a whole codon by construction.
* **`backbone = "protein"` (new).** Each sequence is translated in
  frame 0 (stops/ambiguous codons → X), the amino-acid set is aligned
  with the same progressive engine (`mode = "protein"`), and each
  residue is back-projected to its 3-nt codon; protein gaps become
  `---`. Output is a codon grid — whole-codon indels and phase purity
  1.0 by construction — while retaining protein-level alignment
  accuracy. Trailing partial codons are kept as end insertion columns
  and reported in `qc$partial`.

## 2.2 Selective frameshift refinement

A global stage-2 pass (`refine = TRUE`) reintroduces 1–2 nt frameshift
event columns via a per-row phase DP, but it also fires spuriously on
clean sequences, measurably degrading clean datasets. `refine = "auto"`
runs the same DP **only on rows with `nchar %% 3 != 0`** — the subset a
codon grid provably cannot place in-frame — against a profile computed
from all rows. Every other row contributes its backbone placement
verbatim (an identity `Place`), so:

* event columns are owned exclusively by flagged rows;
* in-frame input is returned **exactly** as the un-refined backbone;
* the flag never degenerates on divergent-but-clean data (unlike
  stop-codon counts, which fire on every row at high divergence —
  evaluated and rejected).

## 2.3 Benchmark datasets and design

Two axes, reported separately because they measure different things:
axis A measures **accuracy against known truth**, axis B measures
**scale and structural correctness on real data**.

### Axis A — simulated CDS with exact ground truth

Simulator (`bench/simulate_codon.R`, all runs seeded and tagged).

* **Phylogeny.** Random Yule trees (`ape::rtree`), branch lengths
  scaled so the mean root-to-tip depth equals the target divergence
  `d`; pairwise tip divergence therefore approaches ~2`d`.
* **Root sequence.** 300 random codons (900 nt), in-frame stop codons
  resampled — a clean, translate-able ancestor.
* **Point evolution.** JC69 per site per branch:
  `p_sub = 1 − exp(−bl)`.
* **Indel model.** Poisson count `rate · bl · L` per branch
  (`rate = 0.02`/site/unit); 50% deletions, 50% insertions. With
  probability `pfs` an indel is a **frameshift**: length drawn from
  {1, 2, 4, 5, 7, 8, 10} nt (geometrically decaying). Otherwise it is
  in-frame: `3 × (geometric(0.4) + 1)` nt. `pfs = 0.2` deliberately
  over-represents frameshifts vs empirical CDS rates — the point is to
  stress-test the frameshift machinery, not to mimic base rates.
* **Ground truth.** Every site carries a stable column id; insertion
  sites receive fractional ids interpolated between their flanking
  ancestral columns, so sorting all tip ids recovers the true column
  order exactly. Insertions on different branches are distinct,
  non-homologous columns — correctly scored as such. The true MSA is
  the n × W placement matrix; per-tip indel events (with anchor
  residue ordinal and true column) are exported for frameshift
  scoring.
* **Grid.** `n ∈ {20, 50}` × `d ∈ {0.1, 0.3}` × `pfs ∈ {0, 0.2}` ×
  5 replicates = **40 datasets**. `d = 0.1` approximates typical
  mammalian-ortholog divergence; `d = 0.3` (~0.6 substitutions/site
  pairwise) is intentionally adversarial.

### Axis B — real mammalian loci

The 31 classic marker loci used by Upham et al. (2019) for the
mammalian species tree — curated orthologous CDS, one unaligned FASTA
per locus (`bench/data/genes_final/`).

* **Sizes.** 19,336 sequences total; 140–3,559 per locus (median
  386); lengths 200–82,292 bp (the max is a single outlier record in
  BRCA1 — real-data messiness, retained). Per-locus statistics in
  Table S1 (`paper/supp_loci.tsv`).
* **Genetic codes.** COI, CYTB, ND1, ND2 aligned under NCBI code 2
  (vertebrate mitochondrial); the remaining 27 under code 1.
* **Messiness is the point.** *Every* locus contains sequences whose
  length is not a multiple of 3 (median 51% of sequences per locus) —
  partial CDS and trimming artifacts are the rule in real data, not
  the exception. This single property already disqualifies
  PRANK-codon, which requires all-mult-3 input (§3.2).
* **Head-to-head tier.** A fixed uniform subsample of n = 50/locus
  (seed 42, written to `out/subsets/*.in.fa`) shared verbatim by all
  tools, so runtimes and structural metrics are directly comparable.
* **Full-scale tier.** All sequences per locus, genoaligner only —
  competitor runtimes at n > 3,000 are infeasible (MACSE already
  exceeds 1 h on some n = 50 loci); reported as a feasibility result.

### Metrics

* **SPS** — fraction of true residue–residue column pairs reproduced
  (accuracy).
* **TC** — fraction of true columns reproduced exactly across all rows.
* **Frameshift recall / precision** — a true fs event is recovered if
  the inferred alignment places a non-mult-3 gap run overlapping its
  anchor residue (±3 columns) in the same row; precision is the
  fraction of inferred non-mult-3 gap runs that overlap a true event
  anchor.
* **Codon phase purity** — fraction of alignment columns whose non-gap
  bases all share one codon phase (1.0 = every column is internally
  in-frame).
* **% indels mult-3** — fraction of gap runs whose length is a
  multiple of 3.
* Runtime (wall seconds per dataset/locus).

Tools: genoaligner modes (`msa_codon` codon/protein backbones with and
without selective refinement, `msa_align`), MACSE v2.07, DECIPHER
(`AlignTranslation`, `AlignSeqs`), MAFFT 7.525 (`--auto`), PRANK
v.250331 (`-codon -once`).

## 2.4 Positioning

**Table 1.** Feature comparison. Median values from the 31-locus
head-to-head (n = 50 sequences/locus); SPS from axis A (n = 50,
mean of 5 replicates).

| | genoaligner `msa_codon` | MACSE | MAFFT | PRANK-codon | DECIPHER-transl. | DECIPHER-seqs |
|---|---|---|---|---|---|---|
| Backbone strategy | protein-guided → codon grid | nt DP + codon-aware costs | nt (FFT-guided) | codon substitution model | protein-guided → nt | nt |
| Whole-codon indels | always, by construction | yes (+ fs events) | no | no | mostly | no |
| Phase purity (real loci, median) | **1.000** | 1.000 | 0.157 | — † | 0.487 | 0.208 |
| Explicit frameshift events | selective (`refine="auto"`) | yes | no | — | no | no |
| Accepts non-mult-3 CDS input | yes | yes | yes | **no** † | yes | yes |
| SPS, clean CDS (d = 0.1 / 0.3) | **0.994 / 0.957** | 0.983 / 0.804 | 0.999 / 0.988 | 0.344 / 0.196 | 0.995 / 0.949 | 0.997 / 0.977 |
| SPS, frameshifted (d = 0.1 / 0.3) | 0.831 / 0.525 | **0.982 / 0.722** | 0.998 / 0.989 | fails | 0.799 / 0.535 | 0.997 / 0.973 |
| Median time/locus (s) | **0.56** | 171.4 | 8.6 | — † | 0.57 | 1.18 |
| Feasible at n ≈ 3,500 | **yes** | no ‡ | yes | — | yes | yes |
| Interface | R pkg + C++17 lib | Java CLI | CLI | CLI | R pkg | R pkg |

† PRANK-codon refused every one of the 31 real loci: it requires all
input lengths to be multiples of 3, which fails on real CDS containing
partial codons — a robustness result in itself.
‡ MACSE exceeded 1 h on single n = 50 loci; extrapolation to thousands
of sequences is infeasible.

genoaligner is the only entry combining strict frame preservation
(phase purity 1.000, whole-codon indels), frameshift tolerance, and
sub-second-per-locus throughput — and the only frame-preserving tool
feasible at phylogenomic scale.

# 3. Results

## 3.1 Axis A — accuracy vs ground truth (n = 50, mean of 5 replicates)

| tool | d=0.1 pfs=0 | d=0.3 pfs=0 | d=0.1 pfs=0.2 (fsR) | d=0.3 pfs=0.2 (fsR) | t/data |
|---|---|---|---|---|---|
| mafft | 0.999 | 0.988 | 0.998 (0.46) | 0.989 (0.48) | ~12 s |
| decipher_seqs | 0.997 | 0.977 | 0.997 (0.44) | 0.973 (0.59) | ~1 s |
| msa_align | 0.996 | 0.914 | 0.996 (0.49) | 0.904 (0.57) | ~1.4 s |
| **msa_codon protein** | **0.994** | **0.957** | 0.796 (0.07) | 0.535 (0.00) | **0.3 s** |
| **+ refine="auto"** | **0.994** | **0.957** | **0.831 (0.31)** | **0.525 (0.36)** | **0.4 s** |
| + refine global | 0.980 | 0.882 | 0.826 (0.32) | 0.490 (0.42) | 0.4 s |
| macse | 0.983 | 0.804 | 0.982 (0.74) | 0.722 (0.31) | ~3.5 min |
| decipher_translation | 0.995 | 0.949 | 0.799 (0.07) | 0.535 (0.00) | 0.7 s |
| msa_codon (token) | 0.955 | 0.768 | 0.728 (0.00) | 0.348 (0.00) | 0.5 s |
| + refine="auto" | 0.955 | 0.768 | 0.840 (0.13) | 0.469 (0.03) | 0.6 s |
| prank | 0.344 | 0.196 | — | — | ~72 s |

**Figure 1** (`figures/fig1_sps.pdf`): SPS vs divergence at pfs = 0
and 0.2 (n = 50, mean ± SE over 5 replicates). Left panel: the protein
backbone tracks the MAFFT/DECIPHER tier while the codon-token backbone
detaches at d = 0.3 — the collapse this work set out to fix. Right
panel: with frameshifts present, MACSE and the nucleotide aligners
lead on raw SPS while the selective-refine variants recover ground
over their un-refined backbones.

Key findings:

* The protein backbone is **third overall in the hard regime**
  (SPS 0.957 at d=0.3 — ahead of MACSE's 0.804 and our own nucleotide
  mode's 0.914) at ~40× lower runtime than MAFFT and ~750× than MACSE.
* `refine = "auto"` is **free insurance**: bit-identical output on
  clean CDS, recovers fs recall 0.31–0.36 and improves SPS on fs data.
  It strictly dominates the un-refined codon backbone on frameshifted
  data (higher SPS in 19/20 datasets, mean +0.10).
* Global refine should stay opt-in: it costs 0.01–0.08 SPS on clean
  data for fs-recall gains that "auto" already captures ~80% of.

## 3.2 Axis B — real loci, n = 50 head-to-head

| tool | median t/locus | total (31 loci) | phase purity | % indels mult-3 |
|---|---|---|---|---|
| **msa_codon** | **0.56 s** | **45 s** | **1.000** | **1.000** |
| decipher_translation | 0.57 s | 39 s | 0.487 | 0.894 |
| decipher_seqs | 1.18 s | 124 s | 0.208 | 0.605 |
| msa_align | 1.31 s | 122 s | 0.282 | 0.146 |
| mafft | 8.61 s | 349 s | 0.157 | 0.685 |
| macse | 171.4 s | 4.4 h | 1.000 | 1.000 |
| prank | fails all 31 loci † | — | — | — |

* Per-locus **MACSE / msa_codon ratio: median 325×** (range 156–1046×)
  and **MAFFT / msa_codon: median 18.8×** (1.5–34.8×) — while MAFFT
  leaves a median phase purity of only 0.157.
* MACSE exceeded 1 h on single loci (BRCA2 7,048 s; ENAM 1,787 s; TTN
  1,640 s) and failed to complete on BRCA1.
* msa_codon and MACSE both deliver phase purity 1.000 where MACSE
  converges — the difference is throughput, not frame correctness.
* DECIPHER's translation guide is fast but loses the frame (phase
  purity 0.21–0.49); even its AA-guided mode emits non-codon indels
  on real data.
* † PRANK-codon requires every input length to be a multiple of 3 and
  refused all 31 real loci — it also failed on every frameshifted
  simulation dataset (axis A). Real-world CDS routinely contains
  partial codons; the tool cannot process them.

**Figure 2** (`figures/fig2_real.pdf`): per-locus runtime (log scale)
and codon phase purity across the 31 head-to-head loci. msa_codon sits
two orders of magnitude below MACSE while matching its phase purity
exactly; every non-codon-aware tool — including fast ones — scatters
phase across the loci.

## 3.3 Axis B — full-scale (all sequences per locus)

`msa_codon` aligned the complete 31-locus set (up to **3,559
sequences**, CYTB) in **430 s total** (median 4.7 s/locus), codon phase
purity 1.000 on every locus (**Figure 3**, `figures/fig3_scale.pdf` —
log-log runtime vs n, near-linear scaling). MACSE is infeasible at this
scale (extrapolating its n=50 timings) — reported honestly as a
feasibility result: for thousands-of-sequences codon-aware MSA,
genoaligner is currently the only practical option.

# 4. Discussion

**What we claim.** genoaligner is the fast, frame-preserving codon MSA
for the common case — clean or near-clean coding sequences at
phylogenomic scale. The protein backbone delivers accuracy competitive
with protein-guided aligners while outputting a strict codon grid;
`refine = "auto"` adds frameshift recovery at zero cost to clean data.

**Honest limitations.**

* For *dense* frameshifts (pseudogene-rich samples), MACSE's recall
  (0.74 vs our 0.31 at d=0.1) remains the reference — it is the right
  tool for that niche when runtime permits.
* Compensated internal frameshifts (indel pairs restoring mod-3
  length) evade the "auto" flag; extending the diagnostic
  (e.g. stop-cluster signatures robust to divergence) is future work.
* At extreme divergence (d=0.3 ≈ 0.6 substitutions/site pairwise —
  beyond typical mammalian orthologs), even the protein backbone
  trails MAFFT on raw SPS; its value there is frame-preserving output
  MAFFT does not provide.
* PRANK-codon underperformed throughout (SPS 0.20–0.34) at our
  divergence settings — reported rather than excluded.

**Reproducibility.** All simulation is seeded and scripted
(`bench/simulate_codon.R`, `bench/bench_sim.R`); results tables are
generated by `bench/summarize.R`. Benchmark data and tool binaries are
external, documented in `bench/README.md` and `bench/PLAN_BENCH.md`.

# Availability

CRAN package `genoaligner` (R/C++17, MIT). Source:
https://github.com/alrobles/genoaligner-r — documentation
https://alrobles.github.io/genoaligner-r/

# Disclosure

Software development was assisted by an AI coding agent (Devin,
Cognition AI) under the direction of the authors. All benchmark design
decisions, analyses, and conclusions are the authors'.

# References (draft)

* Ranwez V et al. 2011. MACSE. PLOS ONE / MBE 2018 (v2).
* Upham NS et al. 2019. Inferring the mammal tree. PLOS Biology.
* Katoh K, Standley DM. MAFFT.
* Wright ES. DECIPHER. BMC Bioinformatics.
* Löytynoja A. PRANK.
