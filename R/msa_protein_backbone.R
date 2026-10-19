# Protein-backbone codon alignment — composition layer for msa_codon().
#
# The codon-token progressive DP loses accuracy at high nucleotide
# divergence; the protein-level engine does not. This path therefore
# builds the alignment backbone at amino-acid level (same engine,
# mode = "protein") and back-projects each residue to its 3 nt codon,
# yielding a codon grid whose indels are whole-codon by construction.
# The optional stage-2 codon_refine pass then repairs reading frame and
# emits MACSE-style 1-2 nt frameshift event columns.
#
# Everything here is R composition over engine entries — the vendored
# C++ is untouched (api.hpp's contract stays the only engine truth).

# Codon index order is the engine's: base-4 with T=0,C=1,A=2,G=3,
# index = b0*16 + b1*4 + b2 (see codon_params in msa_ref.cpp).
.codon_table <- local({
    b <- c("T", "C", "A", "G")
    g <- expand.grid(b, b, b)
    codons <- paste0(g[[3L]], g[[2L]], g[[1L]])   # b1 slowest, b3 fastest
    aa1 <- "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    # NCBI 2 (vertebrate mitochondrial): TGA=W, ATA=M, AGA/AGG=*
    aa2 <- aa1
    substr(aa2, 15L, 15L) <- "W"
    substr(aa2, 35L, 35L) <- "M"
    substr(aa2, 47L, 48L) <- "**"
    list(gc1 = setNames(strsplit(aa1, "")[[1L]], codons),
         gc2 = setNames(strsplit(aa2, "")[[1L]], codons))
})

# Translate one sequence in frame 0 (the MSA must share a single frame:
# a row translated at another offset shifts its codon phase against the
# alignment columns). Returns the protein (stops/ambiguous -> "X"), the
# nt held out of the translation (trailing partial codon), stop count
# and number of ambiguous codons.
codon_tr_frame <- function(s, tab, f = 0L) {
    s <- toupper(s)
    L <- nchar(s)
    nc <- (L - f) %/% 3
    c3 <- if (nc > 0)
        substring(s, f + 3L * (0:(nc - 1L)) + 1L, f + 3L * (1:nc))
    else character(0)
    a <- unname(tab[c3])
    amb <- is.na(a)
    n_stops <- sum(!amb & a == "*")
    a[amb | a == "*"] <- "X"
    list(aa   = paste0(a, collapse = ""),
         raw  = c3,
         n_stops = n_stops,
         n_amb   = sum(amb),
         frame = f,
         lead  = substr(s, 1L, f),
         tail  = substr(s, f + 3L * nc + 1L, L))
}

# Back-project a protein MSA matrix to a nucleotide codon grid.
# Returns a list: core (n x 3W char matrix), lead/tail flank strings.
codon_backproject <- function(pm, trs) {
    n <- nrow(pm); W <- ncol(pm)
    core <- matrix("-", n, 3L * W)
    for (i in seq_len(n)) {
        res <- which(pm[i, ] != "-")
        k <- length(res)
        if (!k) next
        tgt <- rep(3L * res - 2L, each = 3L) + rep(0:2, k)
        src <- rep(3L * seq_len(k) - 2L, each = 3L) + rep(0:2, k)
        core[i, tgt] <- strsplit(paste0(trs[[i]]$raw, collapse = ""), "",
                                 fixed = TRUE)[[1L]][src]
    }
    list(core = core,
         lead = vapply(trs, `[[`, "", "lead"),
         tail = vapply(trs, `[[`, "", "tail"))
}

# The protein-backbone path of msa_codon(): returns a C++-shaped result
# list (status/aligned/frame/stops/partial) for new_msa_result().
msa_codon_protein <- function(seqs, gc, refine, fft) {
    tab <- if (gc == 2L) .codon_table$gc2 else .codon_table$gc1
    trs <- lapply(seqs, codon_tr_frame, tab)
    aa <- vapply(trs, `[[`, "", "aa")
    pa <- msa_align(unname(aa), mode = "protein",
                    iter_refine = fft$iter, fft_band = fft$band,
                    fft_lags = fft$lags, fft_min_rel = fft$min_rel)
    pm <- as.matrix(pa)
    bp <- codon_backproject(pm, trs)
    rows <- apply(bp$core, 1L, paste0, collapse = "")
    if (identical(refine, TRUE)) {
        rows <- unname(codon_refine_cpp(unname(seqs), rows, gc, 1L))
    } else if (identical(refine, "auto")) {
        # selective stage-2: only sequences whose length is not a
        # multiple of 3 may produce frameshift event columns — the
        # protein backbone cannot place them in-frame by construction.
        # Clean-CDS rows keep the backbone verbatim, so "auto" never
        # degrades in-frame input (flagging on stop count would fire
        # on every row at high divergence).
        flag <- nchar(seqs) %% 3L != 0L
        if (any(flag)) {
            sel <- codon_refine_sel_cpp(unname(seqs), rows, flag, gc, 1L)
            rows <- sel$aln
        }
    }
    # flanks are unplaced nt: emit as insertion columns at the ends so no
    # sequence content is dropped (counted in qc$partial).
    maxl <- max(nchar(bp$lead)); maxt <- max(nchar(bp$tail))
    if (maxl > 0L) rows <- paste0(
        str_pad_right(bp$lead, maxl, "-"), rows)
    if (maxt > 0L) rows <- paste0(
        rows, str_pad_right(bp$tail, maxt, "-"))
    list(status = "ok", error = NULL,
         aligned = rows,
         frame = vapply(trs, `[[`, 0L, "frame"),
         stops = vapply(trs, `[[`, 0L, "n_stops"),
         partial = vapply(trs, function(t)
             t$n_amb + (nchar(t$lead) > 0L) + (nchar(t$tail) > 0L),
             numeric(1)))
}

str_pad_right <- function(x, w, pad = "-") {
    paste0(x, strrep(pad, pmax(0L, w - nchar(x))))
}
