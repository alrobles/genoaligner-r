# genoaligner (development version)

# genoaligner 1.0.0

Initial release.

* Edit-distance (Levenshtein / WFA-equivalent) and Smith-Waterman local
  alignment from a self-contained C++17 CPU core; both with score and CIGAR.
* Pipeline-friendly, vectorised `align()` and `align_sw()` over
  data.frame/tibble pair columns; `smax`-bounded resolution that reports
  unresolved pairs rather than returning silently wrong scores.
* Cran-buildable and checkable without a GPU toolchain; `R CMD check` clean.
* Apptainer recipe (`containers/`) to reproduce the build in a pinned cluster
  environment; verified on a KU HPC Q6000 node with `--gres` + `--nv`.