# genoaligner 2.0.0.9000 (development)

* `msa_codon()` gains `backbone = "protein"`: translate each sequence in
  frame 0, build the alignment backbone at amino-acid level with the same
  progressive engine (`msa_align(mode = "protein")`), then back-project
  every residue to its codon — output is a codon grid with whole-codon
  indels by construction. Simulation benchmarks show substantially higher
  sum-of-pairs accuracy than the codon-token backbone at high divergence
  (SPS ≈ 0.95 vs ≈ 0.5–0.7 on n=50, d=0.3 datasets) while preserving
  codon phase purity 1.0. `backbone = "codon"` (default) is unchanged.
  Nucleotides outside complete codons (trailing partial codons) are kept
  as end insertion columns and counted in `qc$partial`; `refine = TRUE`
  runs the stage-2 codon refinement on the back-projected grid.
* `msa_codon(refine = "auto")` — selective stage-2 refinement: only
  sequences whose length is not a multiple of 3 may emit frameshift
  event columns (`codon_refine_sel_cpp` runs the per-row phase DP on
  flagged rows against the full profile; every other row contributes
  its backbone placement verbatim). In-frame input is returned exactly
  as the un-refined backbone, so "auto" cannot degrade clean CDS —
  on simulated frameshift data it improves SPS over the un-refined
  path and recovers frameshift events (mean fs_recall ≈ 0.3–0.5 on
  the protein backbone, ≈ 0.1 on the codon backbone).

# genoaligner 2.0.0

* **Breaking**: `align()` renamed to `align_edit()` — the generic name
  collides with `Biostrings`/`IRanges` exports in mixed sessions. `align()`
  remains as a deprecated alias returning identical results.
* New multiple-sequence-alignment API over the same portable C++17 core
  (vendored host engine, verified bit-exact against the GPU driver in the
  sibling C++ repository):
  - `msa_align(seqs, mode = "dna" | "protein")` — progressive
    profile-profile MSA with a deterministic NJ guide tree.
  - `msa_codon(seqs, gc = 1 | 2, refine = FALSE)` — MACSE-class codon-aware
    MSA: whole-codon indels preserve the reading frame; `refine = TRUE`
    runs a frameshift-refinement pass producing a nucleotide MSA.
* Both return S3 objects of class `genoaligner_msa` (`$aligned` character
  matrix, `$qc` per-sequence QC — frame / stop codons / partial codons in
  codon mode —, `$params`), with `print`/`as.matrix`/`as.character` methods
  and optional `ape::DNAbin` coercion via `as_dnabin()`.
* `msa_align()`/`msa_codon()` gain iterative refinement and FFT anchors
  (MAFFT `FFT-NS-i` class): `iter_refine = k` realigns the two induced
  profiles of every guide-tree edge and keeps a candidate only when the
  sum-of-pairs score strictly improves; `fft_band > 0` restricts each
  realignment DP to a band around Fourier-detected homology anchors
  (faster, exact-rescored, falls back to the full DP on weak signal).
  On 37 real Lepus cytb sequences two rounds improve the SP objective
  620,764 → 628,384 at ~30% less refinement time with `fft_band = 32`.
* Re-vendored engine: genoaligner-devel `devin/msa-iter-refine` @ 2b15a7a
  (PR #13) — see inst/VENDORED.txt.

# genoaligner 1.0.0.9000

* New dataset `lepus_cytb`: 37 real NCBI cytochrome b records for three hare
  species plus an *Oryctolagus* outgroup (mixed lengths, includes two
  mislabelled records used as a QC teaching case).
* New vignettes: `lepus-cytb` (real-data case study: locus QC, common-window
  cropping, barcode gap, `smax` semantics) and `simulated-use-cases`
  (bounded off-target scoring, adapter scan, read dereplication).

# genoaligner 1.0.0

Initial release.

* Edit-distance (Levenshtein / WFA-equivalent) and Smith-Waterman local
  alignment from a self-contained C++17 CPU core; both with score and CIGAR.
* Pipeline-friendly, vectorised `align()` and `align_sw()` over
  data.frame/tibble pair columns; `smax`-bounded resolution that reports
  unresolved pairs rather than returning silently wrong scores.
* CRAN-buildable and checkable without a GPU toolchain; `R CMD check` clean.
* Apptainer recipe (`containers/`) to reproduce the build in a pinned cluster
  environment; verified on a KU HPC Q6000 node with `--gres` + `--nv`.