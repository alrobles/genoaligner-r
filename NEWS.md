# genoaligner 2.0.0.9000 (development)

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