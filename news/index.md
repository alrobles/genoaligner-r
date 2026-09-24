# Changelog

## genoaligner (development version)

- New dataset `lepus_cytb`: 37 real NCBI cytochrome b records for three
  hare species plus an *Oryctolagus* outgroup (mixed lengths, includes
  two mislabelled records used as a QC teaching case).
- New vignettes: `lepus-cytb` (real-data case study: locus QC,
  common-window cropping, barcode gap, `smax` semantics) and
  `simulated-use-cases` (bounded off-target scoring, adapter scan, read
  dereplication).

## genoaligner 1.0.0

CRAN release: 2026-09-24

Initial release.

- Edit-distance (Levenshtein / WFA-equivalent) and Smith-Waterman local
  alignment from a self-contained C++17 CPU core; both with score and
  CIGAR.
- Pipeline-friendly, vectorised
  [`align()`](https://alrobles.github.io/genoaligner-r/reference/align.md)
  and
  [`align_sw()`](https://alrobles.github.io/genoaligner-r/reference/align_sw.md)
  over data.frame/tibble pair columns; `smax`-bounded resolution that
  reports unresolved pairs rather than returning silently wrong scores.
- CRAN-buildable and checkable without a GPU toolchain; `R CMD check`
  clean.
- Apptainer recipe (`containers/`) to reproduce the build in a pinned
  cluster environment; verified on a KU HPC Q6000 node with `--gres` +
  `--nv`.
