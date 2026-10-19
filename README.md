# genoaligner

[![CRAN status](https://www.r-pkg.org/badges/version/genoaligner)](https://CRAN.R-project.org/package=genoaligner)

**Pairwise and multiple sequence alignment for R, from one portable C++17 core.**

`genoaligner` provides edit-distance (Levenshtein / WFA-equivalent) and
Smith-Waterman local pairwise alignment with score and CIGAR reconstruction,
plus progressive multiple sequence alignment in DNA, protein and codon-aware
modes — all from a self-contained C++17 core that needs no GPU toolchain, so
it builds and runs on CRAN and anywhere R runs.

```r
# pairwise, vectorised over pair columns
align_edit(c("AAAACCC", "ACGT"), c("AAAATCC", "ACCT"), smax = 8)
align_sw(c("ACGTACGT", "GGGG"), "ACGTTCGT", scoring = c(2, -3, 5, 2))

# multiple sequence alignment: sequences in -> alignment + QC out
msa_align(c("ACGTACGTAA", "ACGTTCGT", "ACGTACG"))
msa_codon(lepus_cytb$sequence, names = lepus_cytb$accession, gc = 2)
```

## Install

From CRAN:

```r
install.packages("genoaligner")
```

From GitHub (development head):

```r
remotes::install_github("alrobles/genoaligner-r")
```

## What is this?

A **pairwise** aligner *and* a **multiple** aligner over the same core:

- **`align_edit()`** — edit distance with CIGAR, bounded by `smax`. Pairs
  whose true distance exceeds `smax` are reported unresolved (never silently
  wrong).
- **`align_sw()`** — local Smith-Waterman (affine gaps) with CIGAR and the
  aligned span's coordinates.
- **`msa_align()`** — progressive profile–profile MSA with a deterministic
  neighbor-joining guide tree, in `dna` or `protein` mode. Optional
  `iter_refine` rounds realign every guide-tree bipartition, and `fft_band`
  restricts each realignment to a band around Fourier-detected homology
  anchors (MAFFT `FFT-NS-i` class).
- **`msa_codon()`** — codon-aware MSA (MACSE-class): whole-codon indels
  preserve the reading frame by construction; `refine = "auto"` runs a
  selective stage-2 frameshift pass on rows whose internal codon tiling
  leaves held-out blocks (guided re-tokenization evidence — terminal
  partial codons never fire), and `refine = TRUE` runs it on every row.
  NCBI genetic codes 1 (standard) and 2 (vertebrate mitochondrial —
  CYTB/COI/ND).

`msa_*` return a `genoaligner_msa` object carrying the aligned matrix plus
per-sequence QC (reading frame, stop codons, partial codons) — the QC the
library measured travels inside the result. `as_dnabin()` coerces to
`ape::DNAbin`.

> **v2 breaking change:** `align()` was renamed **`align_edit()`** — the
> generic name collides with `Biostrings`/`IRanges` exports in mixed
> sessions. `align()` remains as a deprecated alias returning identical
> results.

The wavefront/GPU backend (ROCm/CUDA) that computes the same results with
acceleration lives in the sibling C++ repository.

## Documentation

See the [site](https://alrobles.github.io/genoaligner-r/): function reference,
the intro vignette, and the codon-aware MSA case study aligning 37 real hare
cytochrome-b records.
