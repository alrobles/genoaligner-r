# genoaligner

[![CRAN
status](https://www.r-pkg.org/badges/version/genoaligner)](https://CRAN.R-project.org/package=genoaligner)

**Pairwise sequence alignment for R, from one portable C++17 core.**

`genoaligner` provides edit-distance (Levenshtein / WFA-equivalent) and
Smith-Waterman local alignment, both with score and CIGAR
reconstruction, from a self-contained C++17 core that needs no GPU
toolchain — so it builds and runs on CRAN and anywhere R runs.

The two entry points are designed to slot straight into a
data.frame/tibble pipeline:

``` r

align(c("AAAACCC", "ACGT"), c("AAAATCC", "ACCT"), smax = 8)
align_sw(c("ACGTACGT", "GGGG"), "ACGTTCGT", scoring = c(2, -3, 5, 2))
```

## Install

From CRAN:

``` r

install.packages("genoaligner")
```

From GitHub (development head):

``` r

remotes::install_github("alrobles/genoaligner-r")
```

## What is this?

A **pairwise** aligner — not a multiple aligner (it does not replace
MAFFT or MACSE) and not a search tool. It aligns two sequences you
already have:

- **[`align()`](https://alrobles.github.io/genoaligner-r/reference/align.md)**
  — edit distance with CIGAR, bounded by `smax`. Pairs whose true
  distance exceeds `smax` are reported unresolved (never silently
  wrong).
- **[`align_sw()`](https://alrobles.github.io/genoaligner-r/reference/align_sw.md)**
  — local Smith-Waterman (affine gaps) with CIGAR and the aligned span’s
  coordinates.

The wavefront/GPU backend (ROCm/CUDA) that computes the same results
with acceleration lives in the sibling C++ repository.

## Documentation

See the package vignette and function reference in the
[site](https://alrobles.github.io/genoaligner-r/).
