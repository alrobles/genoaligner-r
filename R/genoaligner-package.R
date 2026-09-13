#' genoaligner: GPU-portable pairwise sequence alignment
#'
#' Pairwise sequence alignment from one portable C++17 core: edit distance
#' (Levenshtein / WFA-equivalent) and Smith-Waterman local alignment, both with
#' score and CIGAR reconstruction. The core builds and runs anywhere (it is the
#' CPU backend, so it needs no GPU toolchain); the sibling C++ library tracks the
#' GPU (ROCm/CUDA) backend that computes the same results with wavefront
#' acceleration.
#'
#' The two entry points are \code{\link{align}} (edit distance bounded by
#' \code{smax}) and \code{\link{align_sw}} (local Smith-Waterman with affine
#' gaps). Both are vectorised over pairs, so a two-column pipeline
#' (query, reference) maps straight to a data.frame of scores and CIGARs.
#'
#' @name genoaligner-package
#' @aliases genoaligner-package genoaligner
#' @docType package
#' @md
#' @import Rcpp
#' @useDynLib genoaligner, .registration = TRUE
"_PACKAGE"