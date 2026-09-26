# genoaligner — public R API.
#
# Two entry points, exactly like the C++ library:
#   align_edit()  edit distance (Levenshtein / WFA-equivalent), bounded by smax
#   align_sw()    Smith-Waterman local alignment (affine gaps)
# (align() remains as a deprecated alias of align_edit() — the generic name
# collided with Biostrings/IRanges exports in mixed sessions.)
# Both are VECTORISED over the pairs so they slot straight into a
# data.frame/tibble pipeline. This is the "less spaghetti" contract:
# a column of queries and a column of references in, a data.frame of
# scores/CIGARs out, no glue.

#' Align pairs by edit distance
#'
#' Computes the Levenshtein edit distance (substitution = 1, insertion = 1,
#' deletion = 1) between each pair of sequences, with an optional CIGAR
#' reconstruction. This is the same distance the WFA formulation computes;
#' wavefront/GPU speed is a backend concern, the result is identical.
#'
#' A pair is \emph{resolved} only when its true distance is at most \code{smax};
#' otherwise \code{score} is \code{NA} and the pair is reported as
#' \code{resolved = FALSE}. The count of resolved pairs is deliberately visible:
#' a caller that only reads scores can silently miss that most input was
#' abandoned when \code{smax} is set too tightly.
#'
#' @param pattern Character vector; query sequences (recycled to match
#'   \code{text}).
#' @param text Character vector; reference sequences (recycled to match
#'   \code{pattern}). \code{NA} entries produce \code{NA} rows.
#' @param smax Numeric; maximum edit distance searched, in \code{[0, 511]}.
#'   Pairs whose true distance exceeds it are unresolved by design.
#' @param with_cigar Logical; if \code{FALSE} only scores are computed (cheaper,
#'   no traceback) and \code{cigar} is \code{NA}.
#'
#' @return A data.frame with one row per pair and columns:
#'   \code{score} (edit distance, \code{NA} if unresolved), \code{resolved}
#'   (logical), \code{cigar} (over \{M, X, I, D\}, \code{NA} if unresolved or
#'   \code{!with_cigar}), \code{rescore_ok} and \code{wellformed_ok}
#'   (validation flags computed in the library).
#'
#' CIGAR convention (identical for \code{align_edit} and \code{align_sw}):
#' \code{M}/\code{X} consume a text and a pattern base; \code{I} consumes an
#' extra \emph{text} base; \code{D} consumes an extra \emph{pattern} base.
#'
#' @examples
#' align_edit("ACGTACGT", "ACGTTCGT", smax = 8)
#' align_edit(c("AAAACCC", "ACGT"), c("AAAATCC", "ACCT"), smax = 8)
#' @export
align_edit <- function(pattern, text, smax = 64, with_cigar = TRUE) {
    smax <- as.integer(smax)
    if (length(smax) != 1L || is.na(smax) || smax < 0L || smax > 511L)
        stop("smax must be a single integer in [0, 511].")
    stopifnot_character(pattern); stopifnot_character(text)
    n <- max(length(pattern), length(text))
    pattern <- recycle(pattern, n); text <- recycle(text, n)
    score <- integer(n); resolved <- logical(n); cigar <- character(n)
    r_ok <- logical(n); w_ok <- logical(n)
    for (i in seq_len(n)) {
        if (is.na(pattern[i]) || is.na(text[i])) {
            score[i] <- NA_integer_; resolved[i] <- FALSE
            cigar[i] <- NA_character_; r_ok[i] <- FALSE; w_ok[i] <- FALSE
            next
        }
        o <- align_one(pattern[i], text[i], smax, isTRUE(with_cigar))
        score[i]    <- o$score
        resolved[i] <- o$resolved
        cigar[i]    <- if (is.na(o$cigar)) NA_character_ else as.character(o$cigar)
        r_ok[i]     <- o$rescore_ok
        w_ok[i]     <- o$wellformed_ok
    }
    data.frame(score = score, resolved = resolved, cigar = cigar,
               rescore_ok = r_ok, wellformed_ok = w_ok)
}

#' Align pairs by edit distance (deprecated alias)
#'
#' \code{align()} is deprecated in favour of \code{\link{align_edit}()}: the
#' generic name collides with \pkg{Biostrings}/\pkg{IRanges} exports in
#' sessions that load both. It still works and returns identical results;
#' new code should call \code{align_edit()}.
#'
#' @inheritParams align_edit
#' @return See \code{\link{align_edit}}.
#' @keywords internal
#' @export
align <- function(pattern, text, smax = 64, with_cigar = TRUE) {
    .Deprecated("align_edit")
    align_edit(pattern, text, smax = smax, with_cigar = with_cigar)
}

#' Smith-Waterman local alignment
#'
#' Local alignment with affine gap penalties \{match, mismatch, gap_open,
#' gap_extend\}. Penalties are positive and subtracted. A pair with \code{score
#' == 0} has no positive-scoring local alignment (disjoint sequences); all four
#' coordinates are \code{-1} and \code{cigar} is empty.
#'
#' @param text Character vector; row-axis sequences (recycled).
#' @param pattern Character vector; column-axis sequences (recycled).
#' @param scoring Numeric vector of length 4, named or positional:
#'   \code{(match, mismatch, gap_open, gap_extend)}. Recycled if a single
#'   vector is given. The defaults are an example, not a recommendation.
#' @param with_cigar Logical; if \code{FALSE}, only score and end coordinates
#'   are computed (start coordinates and \code{cigar} are \code{NA}/\code{-1}).
#'
#' @return A data.frame with one row per pair and columns: \code{score},
#'   \code{start_i}, \code{start_j}, \code{end_i}, \code{end_j} (0-based
#'   coordinates of the aligned span, \code{-1} when score is 0/unknown),
#'   \code{cigar} (over the span), \code{rescore_ok}, \code{wellformed_ok} and
#'   \code{too_large} (the pair's matrix exceeded the supported size limit).
#'
#' @examples
#' align_sw("TTTACGTGTT", "ACGTGT", scoring = c(2, -3, 5, 2))
#' align_sw(c("ACGTACGT", "GGGG"), c("ACGTTCGT", "ACGT"), c(2, -3, 5, 2))
#' @export
align_sw <- function(text, pattern, scoring = c(2L, -3L, 5L, 2L),
                     with_cigar = TRUE) {
    if (length(scoring) != 4L)
        stop("scoring must be a vector of length 4: (match, mismatch, gap_open, gap_extend).")
    sc <- as.integer(scoring)
    stopifnot_character(text); stopifnot_character(pattern)
    n <- max(length(text), length(pattern))
    text <- recycle(text, n); pattern <- recycle(pattern, n)
    score <- integer(n); start_i <- start_j <- end_i <- end_j <- integer(n)
    cigar <- character(n); r_ok <- w_ok <- t_large <- logical(n)
    for (i in seq_len(n)) {
        if (is.na(text[i]) || is.na(pattern[i])) {
            score[i] <- NA_integer_; start_i[i] <- start_j[i] <- NA_integer_
            end_i[i] <- end_j[i] <- NA_integer_; cigar[i] <- NA_character_
            r_ok[i] <- w_ok[i] <- t_large[i] <- FALSE; next
        }
        o <- align_sw_one(text[i], pattern[i], sc[1], sc[2], sc[3], sc[4],
                          isTRUE(with_cigar))
        score[i] <- o$score; start_i[i] <- o$start_i; start_j[i] <- o$start_j
        end_i[i] <- o$end_i; end_j[i] <- o$end_j
        cigar[i] <- if (is.na(o$cigar)) NA_character_ else as.character(o$cigar)
        r_ok[i] <- o$rescore_ok; w_ok[i] <- o$wellformed_ok
        t_large[i] <- o$too_large
    }
    data.frame(score = score, start_i = start_i, start_j = start_j,
               end_i = end_i, end_j = end_j, cigar = cigar,
               rescore_ok = r_ok, wellformed_ok = w_ok, too_large = t_large)
}

# -- internal helpers ---------------------------------------------------------
stopifnot_character <- function(x) {
    if (!is.character(x)) stop("argument must be character.")
    invisible(x)
}

recycle <- function(x, n) {
    if (length(x) == n) return(x)
    if (length(x) != 1L) stop("arguments must be length 1 or the common length.")
    rep(x, n)
}

seq_len <- function(n) if (n == 0L) integer(0) else seq.int(1L, n)