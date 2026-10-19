#!/usr/bin/env Rscript
# bench_real_variants.R — axis-B re-run for the new msa_codon variants
# on the SAME seeded n=50 subsets (bench/out/subsets/*.in.fa), plus a
# full-scale pass for backbone="protein". Appends to
# bench/out/bench_results.tsv / bench_results_fullscale.tsv.
#
# Tools: msa_codon_protein (backbone="protein"),
#        msa_codon_protein_auto (backbone="protein", refine="auto"),
#        msa_codon_auto (codon backbone, refine="auto")

suppressPackageStartupMessages(library(genoaligner))

BENCH <- normalizePath("bench"); OUT <- file.path(BENCH, "out")
SUB   <- file.path(OUT, "subsets")
DATA  <- file.path(BENCH, "data", "genes_final")
TSV   <- file.path(OUT, "bench_results.tsv")
TSVFS <- file.path(OUT, "bench_results_fullscale.tsv")
MITO  <- c("COI", "CYTB", "ND1", "ND2")

logm <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
                          sprintf(...), "\n")

read_fasta <- function(path) {
    l <- readLines(path); h <- grep("^>", l)
    ends <- c(h[-1] - 1L, length(l))
    s <- vapply(seq_along(h), function(i)
        paste(l[(h[i] + 1L):ends[i]], collapse = ""), character(1))
    names(s) <- sub("^>\\s*", "", l[h]); toupper(s)
}

indel_stats <- function(aln_rows) {
    lens <- unlist(lapply(aln_rows, function(s) {
        g <- gregexpr("-+", s)[[1]]
        if (g[1] < 0) integer(0) else attr(g, "match.length")
    }))
    if (!length(lens)) return(c(runs = 0, pct_mult3 = NA_real_))
    c(runs = length(lens), pct_mult3 = mean(lens %% 3L == 0L))
}

phase_purity <- function(aln_rows) {
    m <- do.call(rbind, strsplit(aln_rows, "", fixed = TRUE))
    phase <- t(vapply(seq_len(nrow(m)), function(i) {
        isbase <- m[i, ] != "-"
        (cumsum(isbase) - 1L) %% 3L
    }, integer(ncol(m))))
    col_pure <- vapply(seq_len(ncol(m)), function(j) {
        ph <- phase[m[, j] != "-", j]
        length(unique(ph)) <= 1L
    }, logical(1))
    mean(col_pure)
}

record <- function(tsv, ...) {
    row <- data.frame(..., stringsAsFactors = FALSE)
    write.table(row, tsv, sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = !file.exists(tsv), append = file.exists(tsv))
}

run_variant <- function(tool, seqs, gc) {
    t <- system.time({
        a <- switch(tool,
            msa_codon_protein      = msa_codon(seqs, gc = gc,
                                             backbone = "protein"),
            msa_codon_protein_auto = msa_codon(seqs, gc = gc,
                                               backbone = "protein",
                                               refine = "auto"),
            msa_codon_auto         = msa_codon(seqs, gc = gc,
                                               refine = "auto"))
    })
    list(rows = apply(as.matrix(a), 1L, paste0, collapse = ""),
         t = t[["elapsed"]])
}

tools <- c("msa_codon_protein", "msa_codon_protein_auto", "msa_codon_auto")

## ---- tier 1: n=50 head-to-head subsets --------------------------------------
done <- if (file.exists(TSV)) read.delim(TSV) else data.frame()
for (f in sort(list.files(SUB, pattern = "\\.in\\.fa$"))) {
    locus <- sub("\\.in\\.fa$", "", f)
    seqs <- read_fasta(file.path(SUB, f))
    gc <- if (locus %in% MITO) 2L else 1L
    for (tool in tools) {
        if (nrow(done) && any(done$locus == locus & done$tool == tool)) next
        r <- tryCatch(run_variant(tool, seqs, gc),
                      error = function(e) list(rows = character(0), t = NA))
        if (!length(r$rows)) {
            record(TSV, locus = locus, n_seqs = length(seqs), tool = tool,
                   t_sec = NA, aln_cols = NA, indel_runs = NA,
                   pct_indel_mult3 = NA, codon_phase_purity = NA,
                   status = "error")
            logm("%s %s ERROR", locus, tool); next
        }
        is <- indel_stats(r$rows); pp <- phase_purity(r$rows)
        record(TSV, locus = locus, n_seqs = length(seqs), tool = tool,
               t_sec = round(r$t, 2), aln_cols = nchar(r$rows[1]),
               indel_runs = unname(is["runs"]),
               pct_indel_mult3 = round(unname(is["pct_mult3"]), 4),
               codon_phase_purity = round(pp, 4), status = "ok")
        logm("%s %-24s %6.1fs cols=%d m3=%.2f phase=%.3f",
             locus, tool, r$t, nchar(r$rows[1]), is["pct_mult3"], pp)
    }
}

## ---- tier 2: full scale, protein backbone ------------------------------------
tools_fs <- c("msa_codon_protein", "msa_codon_protein_auto")
done <- if (file.exists(TSVFS)) read.delim(TSVFS) else data.frame()
for (f in sort(list.files(DATA, pattern = "\\.fasta$"))) {
    locus <- sub("\\.fasta$", "", f)
    seqs <- read_fasta(file.path(DATA, f))
    gc <- if (locus %in% MITO) 2L else 1L
    for (tool in tools_fs) {
        if (nrow(done) && any(done$locus == locus & done$tool == tool)) next
        r <- tryCatch(run_variant(tool, seqs, gc),
                      error = function(e) list(rows = character(0), t = NA))
        if (!length(r$rows)) {
            record(TSVFS, locus = locus, n_seqs = length(seqs), tool = tool,
                   t_sec = NA, aln_cols = NA, indel_runs = NA,
                   pct_indel_mult3 = NA, codon_phase_purity = NA,
                   status = "error")
            logm("FS %s %s ERROR", locus, tool); next
        }
        is <- indel_stats(r$rows); pp <- phase_purity(r$rows)
        record(TSVFS, locus = locus, n_seqs = length(seqs), tool = tool,
               t_sec = round(r$t, 2), aln_cols = nchar(r$rows[1]),
               indel_runs = unname(is["runs"]),
               pct_indel_mult3 = round(unname(is["pct_mult3"]), 4),
               codon_phase_purity = round(pp, 4), status = "ok")
        logm("FS %s %-24s %6.1fs cols=%d m3=%.2f phase=%.3f n=%d",
             locus, tool, r$t, nchar(r$rows[1]), is["pct_mult3"], pp,
             length(seqs))
    }
}
logm("done")
