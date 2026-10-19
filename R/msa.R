# genoaligner — multiple sequence alignment (MSA) R API.
#
# Two entry points over the same C++ engine (progressive profile-profile
# alignment with an NJ guide tree; the vendored host engine is verified
# bit-exact against the GPU driver in the C++ repository):
#   msa_align()  DNA / protein MSA
#   msa_codon()  codon-aware MSA (MACSE-class): whole-codon indels by default,
#                optional stage-2 frameshift refinement
#
# Both return an S3 object of class "genoaligner_msa"; the QC the library
# measured during encoding travels inside it ($qc), same spirit as the
# resolved/flags columns of the pairwise API.

#' Multiple sequence alignment (DNA or protein)
#'
#' Progressive MSA: fragment-corrected k-mer distances, a deterministic
#' neighbor-joining guide tree, post-order profile-vs-profile semiglobal
#' alignment (free end gaps, position-specific gap penalties). The pipeline
#' is deterministic: identical input always produces identical output.
#'
#' @param seqs Character vector; unaligned sequences (no \code{"-"}).
#'   IUPAC nucleotide letters in \code{dna} mode; amino-acid letters in
#'   \code{protein} mode.
#' @param names Optional character vector of sequence names, same length as
#'   \code{seqs}; defaults to \code{names(seqs)} or generated ids. Used for
#'   row names of the aligned matrix and FASTA output.
#' @param mode \code{"dna"} (default) or \code{"protein"}.
#' @param guide Guide-tree construction; only \code{"nj"} (deterministic
#'   neighbor joining) exists today.
#' @param guide_tree Optional \code{ape::phylo} guide tree overriding
#'   \code{guide}: the progressive merges follow its topology (branch
#'   lengths are ignored). Must be fully bifurcating with
#'   \code{tip.label} matching the sequence names. Diagnostic arm for
#'   benchmarking — e.g. pass the simulation's true tree to separate
#'   guide-tree error from scoring error.
#' @param iter_refine Integer \eqn{\ge 0}; iterative-refinement rounds over
#'   the guide-tree bipartitions (MAFFT \code{FFT-NS-i} class): every tree
#'   edge splits the alignment into two induced sub-profiles, the pair is
#'   realigned, and the candidate is kept only if the sum-of-pairs score
#'   strictly improves. \code{0} (default) keeps the plain progressive
#'   alignment. Two rounds is a good first value; the pass converges and a
#'   further round reports no accepted changes.
#' @param fft_band Integer \eqn{\ge 0}; FFT anchor band half-width in
#'   columns. When \eqn{> 0} (and \code{iter_refine > 0}) the engine
#'   detects homology anchors by Fourier cross-correlation of the two
#'   profile property signals and runs the realignment DP inside a band
#'   around them — the literal FFT of MAFFT \code{FFT-NS-i}. Weak anchor
#'   evidence falls back to the full DP, and the strict-improvement gate
#'   applies either way, so this is a speed/robustness knob, not a
#'   quality risk. \code{0} (default) uses the full matrix.
#' @param fft_lags Integer \eqn{\ge 1}; number of anchor diagonals whose
#'   bands are unioned (default \code{4}). Only used when
#'   \code{fft_band > 0}.
#' @param fft_min_rel Numeric \eqn{\ge 0}; minimum cosine similarity of
#'   the best anchor before the band is trusted (default \code{0.10}).
#'   Higher values make the engine fall back to the full DP more often.
#'
#' @return An object of class \code{genoaligner_msa} with components
#'   \code{aligned} (character matrix, one row per input sequence, one column
#'   per alignment position), \code{qc} (per-sequence data.frame; all
#'   \code{NA} outside codon mode) and \code{params} (echo of the effective
#'   request). See \code{\link{msa_codon}} for the codon-aware mode.
#'
#' @examples
#' msa_align(c("ACGTACGTAA", "ACGTTCGT", "ACGTACG"))
#' msa_align(c(insulin = "MALWMRLLPLL", insulin2 = "MALWTRLLPLL"),
#'           mode = "protein")
#' @export
msa_align <- function(seqs, names = NULL, mode = c("dna", "protein"),
                      guide = "nj", guide_tree = NULL,
                      iter_refine = 0L, fft_band = 0L,
                      fft_lags = 4L, fft_min_rel = 0.10) {
    mode <- match.arg(mode)
    if (!identical(guide, "nj"))
        stop("only guide = \"nj\" is supported.")
    fft <- check_fft_args(iter_refine, fft_band, fft_lags, fft_min_rel)
    seqs <- check_msa_seqs(seqs)
    names <- check_msa_names(names, seqs)
    o <- if (!is.null(guide_tree)) {
        mg <- phylo_to_merges(guide_tree, names)
        msa_tree_cpp(unname(seqs), mode, gc_def = 1L,
                     codon_refine = 0L, local_frame = FALSE,
                     iter_refine = fft$iter, fft_band = fft$band,
                     fft_lags = fft$lags, fft_min_rel = fft$min_rel,
                     mg$mleft, mg$mright)
    } else {
        msa_run_cpp(unname(seqs), mode, gc_def = 1L,
                    codon_refine = 0L, local_frame = FALSE,
                    iter_refine = fft$iter, fft_band = fft$band,
                    fft_lags = fft$lags, fft_min_rel = fft$min_rel)
    }
    new_msa_result(o, seqs, names,
                   params = c(list(mode = mode,
                                   guide = if (is.null(guide_tree))
                                       guide else "provided"),
                              fft$params))
}

#' Codon-aware multiple sequence alignment
#'
#' Codon-aware MSA of coding sequences: each sequence is tokenized to codons
#' (NCBI genetic-code table \code{gc}) and aligned over a codon
#' substitution matrix (BLOSUM62 on the translated amino acids plus a
#' per-identical-nucleotide bonus and a stop-codon penalty). With the default
#' \code{refine = FALSE} every indel is a whole codon, so every aligned row
#' stays a multiple of 3 and the reading frame is preserved by construction.
#'
#' With \code{refine = TRUE} a MACSE-style stage-2 pass reintroduces
#' internal 1-2 nt frameshift events as narrow columns, so the output is a
#' nucleotide MSA no longer guaranteed to be a multiple of 3 — that is the
#' intended semantics for sequences that truly carry frameshifts.
#'
#' @inheritParams msa_align
#' @param gc Integer; NCBI genetic code: \code{1L} (standard, default) or
#'   \code{2L} (vertebrate mitochondrial — the right table for CYTB/COI/ND).
#' @param refine Frameshift refinement. \code{FALSE} (default): no stage-2
#'   pass; \code{TRUE}: run it on every row (helps sequences with true
#'   frameshifts but costs accuracy on clean CDS); \code{"auto"}:
#'   selective — only sequences whose length is not a multiple of 3 may
#'   produce frameshift event columns, so in-frame input is returned
#'   exactly as the un-refined backbone.
#' @param backbone How the alignment backbone is built. \code{"codon"}
#'   (default) tokenizes codons and aligns over the codon substitution
#'   matrix. \code{"protein"} translates each sequence in frame 0,
#'   aligns at amino-acid level (\code{\link{msa_align}} with
#'   \code{mode = "protein"}), then back-projects to a codon grid — more
#'   accurate at high nucleotide divergence, with the same whole-codon
#'   indel guarantee. Nucleotides outside the chosen frame (leading or
#'   trailing partial codons) are kept as insertion columns at the
#'   alignment ends and counted in \code{qc$partial}.
#'
#' @return A \code{genoaligner_msa} object. In codon mode \code{qc} carries
#'   measured values per input sequence: \code{frame} (chosen reading frame,
#'   0-2), \code{n_stops} (stop codons in the tokenized stream) and
#'   \code{partial} (partial/ambiguous codons and held-out frameshift
#'   blocks). A clean coding sequence reports \code{n_stops = 0}.
#'
#' @examples
#' msa_codon(c("ATGATAATCACC", "ATGATTATCACCTGA"))
#' msa_codon(c("ATGATAATCACC", "ATGATTATCACCTGA"), gc = 2)
#' @export
msa_codon <- function(seqs, names = NULL, gc = 1L, refine = FALSE,
                      backbone = c("codon", "protein"),
                      guide_tree = NULL,
                      iter_refine = 0L, fft_band = 0L,
                      fft_lags = 4L, fft_min_rel = 0.10) {
    gc <- as.integer(gc)
    if (length(gc) != 1L || is.na(gc) || !gc %in% c(1L, 2L))
        stop("gc must be a single integer: 1 (standard) or 2 (vertebrate mitochondrial).")
    backbone <- match.arg(backbone)
    refine <- if (is.character(refine)) match.arg(refine, "auto")
              else isTRUE(refine)
    fft <- check_fft_args(iter_refine, fft_band, fft_lags, fft_min_rel)
    seqs <- check_msa_seqs(seqs)
    names <- check_msa_names(names, seqs)
    mg <- if (!is.null(guide_tree)) phylo_to_merges(guide_tree, names)
    o <- if (backbone == "protein") {
        msa_codon_protein(seqs, gc, refine, fft, guide_tree, names)
    } else if (!is.null(mg)) {
        # Supplied guide tree: same dispatch as msa_run_cpp but the
        # progressive stage follows the given topology. "auto" still
        # needs the strict-tokenization grid + selective stage-2.
        oo <- msa_tree_cpp(unname(seqs), "codon", gc_def = gc,
                           codon_refine = if (identical(refine, "auto"))
                               0L else as.integer(refine),
                           local_frame = identical(refine, TRUE),
                           iter_refine = fft$iter, fft_band = fft$band,
                           fft_lags = fft$lags, fft_min_rel = fft$min_rel,
                           mg$mleft, mg$mright)
        if (oo$status == "ok" && identical(refine, "auto") &&
            any(nchar(seqs) %% 3L != 0L)) {
            sel <- codon_refine_sel_cpp(unname(seqs), oo$aligned,
                                        nchar(seqs) %% 3L != 0L, gc, 1L)
            oo$aligned <- sel$aln
        }
        oo
    } else if (identical(refine, "auto")) {
        # Strict frame-0 tokenization gives a pure triplet grid
        # (local_frame cannot be used here: the engine forces an internal
        # refine round on it, and its partial fills break re-tokenization).
        # Then selective stage-2: only rows whose length is not a
        # multiple of 3 may emit frameshift event columns.
        oo <- msa_run_cpp(unname(seqs), "codon", gc_def = gc,
                          codon_refine = 0L, local_frame = FALSE,
                          iter_refine = fft$iter, fft_band = fft$band,
                          fft_lags = fft$lags, fft_min_rel = fft$min_rel)
        if (oo$status == "ok" && any(nchar(seqs) %% 3L != 0L)) {
            sel <- codon_refine_sel_cpp(unname(seqs), oo$aligned,
                                        nchar(seqs) %% 3L != 0L, gc, 1L)
            oo$aligned <- sel$aln
        }
        oo
    } else {
        msa_run_cpp(unname(seqs), "codon", gc_def = gc,
                    codon_refine = if (refine) 1L else 0L,
                    local_frame = refine,
                    iter_refine = fft$iter, fft_band = fft$band,
                    fft_lags = fft$lags, fft_min_rel = fft$min_rel)
    }
    new_msa_result(o, seqs, names,
                   params = c(list(mode = "codon", gc = gc, refine = refine,
                                   backbone = backbone), fft$params))
}

#' Sum-of-pairs score of a finished alignment
#'
#' The objective iterative refinement maximizes internally
#' (\code{iter_refine > 0}): every unordered row pair contributes its
#' induced pairwise alignment, scored through the mode's substitution
#' table minus gap-run costs. Diagnostic use — e.g. checking whether the
#' model's own SP ranks backbone variants the way an external accuracy
#' metric does.
#'
#' For \code{mode = "codon"} the input must be a pure codon grid (row
#' width a multiple of 3, every column a triplet). To score a refined
#' alignment that contains frameshift event columns, subset the codon
#' columns first (the \code{cofs} offsets returned by the selective
#' refinement layer).
#'
#' @param x A \code{genoaligner_msa} object, a character matrix, or a
#'   character vector of equal-length aligned rows.
#' @param mode \code{"dna"} (default), \code{"protein"} or \code{"codon"}.
#' @param gc Integer; NCBI genetic code for \code{mode = "codon"}.
#' @return Numeric scalar; higher is better.
#' @examples
#' a <- msa_align(c("ACGTACGT", "ACGTTCGT"))
#' msa_sp_score(a)
#' @export
msa_sp_score <- function(x, mode = c("dna", "protein", "codon"), gc = 1L) {
    mode <- match.arg(mode)
    rows <- if (inherits(x, "genoaligner_msa"))
        apply(x$aligned, 1L, paste0, collapse = "")
    else if (is.character(x) || is.matrix(x))
        apply(as.matrix(x), 1L, paste0, collapse = "")
    else stop("x must be a genoaligner_msa, matrix, or character vector.")
    gc <- as.integer(gc)
    if (mode == "codon" && (length(gc) != 1L || !gc %in% c(1L, 2L)))
        stop("gc must be 1 or 2 for mode = \"codon\".")
    msa_sp_score_cpp(unname(rows), mode, gc)
}

#' Coerce a genoaligner_msa alignment to ape::DNAbin
#'
#' Requires the \pkg{ape} package (in \code{Suggests}; used here only).
#'
#' @param x A \code{genoaligner_msa} object.
#' @return An \code{ape::DNAbin} matrix.
#' @examples
#' if (requireNamespace("ape", quietly = TRUE)) {
#'     as_dnabin(msa_align(c("ACGTACGT", "ACGTTCGT")))
#' }
#' @export
as_dnabin <- function(x) {
    if (!inherits(x, "genoaligner_msa"))
        stop("x must be a genoaligner_msa object.")
    if (!requireNamespace("ape", quietly = TRUE))
        stop("package 'ape' is needed for DNAbin coercion.")
    ape::as.DNAbin(as.matrix(x))
}

# -- S3 methods --------------------------------------------------------------

#' @export
print.genoaligner_msa <- function(x, ...) {
    p <- x$params
    cat(sprintf("genoaligner MSA: %d sequences x %d columns (mode=%s%s)\n",
                nrow(x$aligned), ncol(x$aligned), p$mode,
                if (p$mode == "codon")
                    sprintf(", gc=%d, refine=%s", p$gc, p$refine)
                else ""))
    if (!all(is.na(x$qc$n_stops))) {
        ns <- x$qc$n_stops
        cat(sprintf("  qc: %d/%d sequences without stop codons",
                    sum(ns == 0L, na.rm = TRUE), sum(!is.na(ns))))
        np <- x$qc$partial
        if (any(np > 0L, na.rm = TRUE))
            cat(sprintf("; %d with partial/frameshift blocks",
                        sum(np > 0L, na.rm = TRUE)))
        cat("\n")
    }
    invisible(x)
}

#' @export
as.matrix.genoaligner_msa <- function(x, ...) {
    x$aligned
}

#' @export
as.character.genoaligner_msa <- function(x, ...) {
    n <- nrow(x$aligned)
    nm <- rownames(x$aligned)
    rows <- apply(x$aligned, 1L, paste0, collapse = "")
    vapply(seq_len(n), function(i)
        paste0(">", nm[i], "\n", rows[i]), character(1))
}

# -- internals ---------------------------------------------------------------

check_msa_seqs <- function(seqs) {
    if (is.null(seqs) || !is.character(seqs) || length(seqs) == 0L)
        stop("seqs must be a non-empty character vector.")
    if (anyNA(seqs))
        stop("seqs must not contain NA.")
    seqs
}

# Convert an ape::phylo to the engine's merge list: leaves are 0-based
# sequence indices (tip.label matched against names), internal nodes
# take ids n + postorder rank so every child's id is smaller than its
# parent's -- the invariant msa_iter_refine/tree_levels rely on.
# Returned mleft[k]/mright[k] are the two children of internal node
# n+k in that order; the last merge is the root join.
phylo_to_merges <- function(phy, names) {
    if (!inherits(phy, "phylo"))
        stop("guide_tree must be an ape::phylo object.")
    n <- length(names)
    if (is.null(phy$edge) || length(phy$tip.label) != n)
        stop("guide_tree must have exactly one tip per sequence.")
    leaf <- match(phy$tip.label, names)
    if (anyNA(leaf) || anyDuplicated(phy$tip.label))
        stop("guide_tree tip labels must match the sequence names exactly.")
    nnode <- phy$Nnode
    if (is.null(nnode) || nnode != n - 1L)
        stop("guide_tree must be fully bifurcating.")
    ch <- vector("list", n + nnode)
    for (e in seq_len(nrow(phy$edge)))
        ch[[phy$edge[e, 1L]]] <- c(ch[[phy$edge[e, 1L]]], phy$edge[e, 2L])
    if (any(lengths(ch[n + seq_len(nnode)]) != 2L))
        stop("guide_tree must be fully bifurcating.")
    root <- setdiff(n + seq_len(nnode), phy$edge[, 2L])
    if (length(root) != 1L)
        stop("guide_tree has no unique root.")
    # DFS preorder over ape ids, reversed: descendants precede ancestors.
    post <- integer(0); st <- root
    while (length(st)) {
        u <- st[length(st)]; st <- st[-length(st)]
        post <- c(post, u)
        if (u > n) st <- c(st, ch[[u]])
    }
    internal <- rev(post[post > n])
    id <- integer(n + nnode)
    id[seq_len(n)] <- leaf - 1L
    id[internal] <- n + seq_along(internal) - 1L
    kids <- ch[internal]
    list(mleft  = id[vapply(kids, `[`, 0L, 1L)],
         mright = id[vapply(kids, `[`, 0L, 2L)])
}

# Refine/FFT request validation at the R boundary — same refusal rules
# as the C++ request gate (the library validates again; this layer keeps
# the error phrasing R-idiomatic).
check_fft_args <- function(iter_refine, fft_band, fft_lags, fft_min_rel) {
    iter <- suppressWarnings(as.integer(iter_refine))
    band <- suppressWarnings(as.integer(fft_band))
    lags <- suppressWarnings(as.integer(fft_lags))
    mr   <- suppressWarnings(as.numeric(fft_min_rel))
    if (length(iter) != 1L || is.na(iter) || iter < 0L)
        stop("iter_refine must be a single integer >= 0.")
    if (length(band) != 1L || is.na(band) || band < 0L)
        stop("fft_band must be a single integer >= 0.")
    if (length(lags) != 1L || is.na(lags) || lags < 1L)
        stop("fft_lags must be a single integer >= 1.")
    if (length(mr) != 1L || !is.finite(mr) || mr < 0)
        stop("fft_min_rel must be a single finite number >= 0.")
    if (band > 0L && iter < 1L)
        stop("fft_band only applies during iterative refinement: ",
             "set iter_refine >= 1 (or fft_band = 0).")
    list(iter = iter, band = band, lags = lags, min_rel = mr,
         params = list(iter_refine = iter, fft_band = band,
                       fft_lags = lags, fft_min_rel = mr))
}

check_msa_names <- function(names, seqs) {
    if (is.null(names)) names <- names(seqs)
    if (is.null(names)) names <- paste0("seq", seq_len(length(seqs)))
    if (!is.character(names) || length(names) != length(seqs))
        stop("names must be a character vector the same length as seqs.")
    names
}

new_msa_result <- function(o, seqs, names, params) {
    if (o$status != "ok") {
        if (o$status == "invalid_argument")
            stop(sub("^msa_align: ", "", o$error), call. = FALSE)
        stop("msa_align failed: ", o$error, call. = FALSE)
    }
    aln <- strsplit(o$aligned, "", fixed = TRUE)
    width <- if (length(aln)) nchar(o$aligned[1]) else 0L
    m <- matrix(unlist(aln), nrow = length(aln), ncol = width,
                byrow = TRUE, dimnames = list(names, NULL))
    structure(
        list(aligned = m,
             qc = data.frame(frame = o$frame, n_stops = o$stops,
                             partial = o$partial),
             params = c(params, list(n = length(seqs), width = width))),
        class = "genoaligner_msa")
}
