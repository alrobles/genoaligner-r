#!/usr/bin/env Rscript
# Phase-5 benchmark: msa_codon / msa_align vs MACSE v2.07 vs DECIPHER
# on the 31 phylogenyAI "classic loci" (real mammalian coding genes).
#
# Usage: Rscript bench/bench_codon_msa.R [--n=N] [--loci=REGEX]
#   --n=N        sequences sampled per locus (0 = all; default 200)
#   --loci=RE    only loci whose name matches RE (default all 31)
#
# Outputs bench/out/bench_results.tsv — one row per (locus, tool):
#   locus, n_seqs, tool, t_sec, aln_cols, indel_runs, pct_indel_mult3,
#   codon_phase_purity, status
#
# Frame-conservation metrics (tool-agnostic):
#   * pct_indel_mult3 — % of per-sequence gap *runs* whose length is a
#     multiple of 3 (whole-codon indels preserve the reading frame).
#   * codon_phase_purity — fraction of alignment columns in which every
#     non-gap base occupies the same codon phase (pos-1/2/3) within its
#     own reading frame. Nucleotide aligners mix phases whenever an
#     indel length is not a multiple of 3.
#
# External tools are invoked, never bundled: MACSE runs from
# bench/tools/macse_v2.07.jar; DECIPHER comes from Bioconductor.

suppressPackageStartupMessages({
    library(genoaligner)
    library(Biostrings)
    library(DECIPHER)
})

args <- commandArgs(TRUE)
opt <- function(key, default) {
    m <- grep(paste0("^--", key, "="), args, value = TRUE)
    if (length(m)) sub(paste0("^--", key, "="), "", m[1]) else default
}
n_per_locus <- as.integer(opt("n", "200"))
loci_re     <- opt("loci", ".")
tools_keep  <- strsplit(opt("tools",
    "msa_codon,msa_align,macse,decipher_translation,decipher_seqs"),
    ",")[[1]]

BENCH <- normalizePath("bench")
DATA  <- file.path(BENCH, "data", "genes_final")
MACSE <- file.path(BENCH, "tools", "macse_v2.07.jar")
OUT   <- file.path(BENCH, "out"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SUB   <- file.path(OUT, "subsets"); dir.create(SUB, showWarnings = FALSE, recursive = TRUE)
TSV   <- file.path(OUT, paste0("bench_results", opt("tag", ""), ".tsv"))
LOG   <- file.path(OUT, "bench.log")

MITO <- c("COI", "CYTB", "ND1", "ND2")  # NCBI code 2 (vertebrate mitochondrial)

logm <- function(...) {
    msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
    cat(msg, "\n"); cat(msg, "\n", file = LOG, append = TRUE)
}

read_fasta <- function(path) {
    ln <- readLines(path, warn = FALSE)
    hdr <- grepl("^>", ln)
    nm <- sub("^>\\s*", "", ln[hdr])
    seqs <- vapply(split(ln[!hdr], cumsum(hdr)[!hdr]),
                   paste0, character(1), collapse = "")
    names(seqs) <- sub("\\s.*$", "", nm)
    toupper(seqs)
}

write_fasta <- function(v, path) {
    con <- file(path, "wt")
    for (i in seq_along(v))
        writeLines(c(paste0(">", names(v)[i]), v[i]), con)
    close(con)
    invisible(path)
}

## ---- metrics ---------------------------------------------------------------

# per-row gap runs -> lengths
indel_stats <- function(aln_rows) {
    lens <- unlist(lapply(aln_rows, function(s) {
        g <- gregexpr("-+", s)[[1]]
        if (g[1] < 0) integer(0) else attr(g, "match.length")
    }))
    if (!length(lens)) return(c(runs = 0, pct_mult3 = NA_real_))
    c(runs = length(lens), pct_mult3 = mean(lens %% 3L == 0L))
}

# fraction of columns whose non-gap bases share one codon phase
phase_purity <- function(aln_rows) {
    m <- do.call(rbind, strsplit(aln_rows, "", fixed = TRUE))
    # phase of every cell = (number of non-gap bases before it in its row) %% 3
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

## ---- tool runners ----------------------------------------------------------

run_genoaligner <- function(seqs, gc, codon = TRUE) {
    t <- system.time({
        a <- if (codon) msa_codon(seqs, names = names(seqs), gc = gc)
             else       msa_align(seqs, names = names(seqs))
    })
    list(rows = apply(as.matrix(a), 1L, paste0, collapse = ""), t = t[["elapsed"]])
}

run_decipher <- function(seqs, gc, translation = FALSE) {
    ss <- DNAStringSet(seqs)
    t <- system.time({
        a <- if (translation)
            AlignTranslation(ss, geneticCode = getGeneticCode(as.character(gc)),
                             verbose = FALSE)
        else
            AlignSeqs(ss, verbose = FALSE)
    })
    list(rows = toupper(as.character(a)), t = t[["elapsed"]])
}

run_macse <- function(seqs, gc, tag) {
    ifile <- file.path(SUB, paste0(tag, ".in.fa"))
    ofile <- file.path(SUB, paste0(tag, ".macse.fa"))
    write_fasta(seqs, ifile)
    t <- system.time({
        system2("java", c("-jar", shQuote(MACSE), "-prog", "alignSequences",
                          "-seq", shQuote(ifile), "-gc_def", gc,
                          "-out_NT", shQuote(ofile)),
                stdout = FALSE, stderr = FALSE)
    })
    if (!file.exists(ofile)) return(list(rows = character(0), t = t[["elapsed"]],
                                        err = TRUE))
    list(rows = read_fasta(ofile), t = t[["elapsed"]])
}

## ---- driver ----------------------------------------------------------------

tools <- intersect(c("msa_codon", "msa_align", "macse",
                     "decipher_translation", "decipher_seqs"), tools_keep)

record <- function(...) {
    row <- data.frame(..., stringsAsFactors = FALSE)
    write.table(row, TSV, sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = !file.exists(TSV), append = file.exists(TSV))
}

files <- list.files(DATA, pattern = "\\.fasta$")
loci <- sub("\\.fasta$", "", files)
loci <- loci[grepl(loci_re, loci)]
logm("benchmark start: %d loci, n=%s", length(loci),
     if (n_per_locus > 0) as.character(n_per_locus) else "all")

for (locus in loci) {
    seqs <- read_fasta(file.path(DATA, paste0(locus, ".fasta")))
    if (n_per_locus > 0 && length(seqs) > n_per_locus) {
        set.seed(42)
        seqs <- seqs[sample.int(length(seqs), n_per_locus)]
    }
    gc <- if (locus %in% MITO) 2L else 1L
    logm("%s: %d seqs, gc=%d", locus, length(seqs), gc)

    for (tool in tools) {
        r <- tryCatch(switch(tool,
            msa_codon             = run_genoaligner(seqs, gc, TRUE),
            msa_align             = run_genoaligner(seqs, gc, FALSE),
            macse                 = run_macse(seqs, gc, locus),
            decipher_translation  = run_decipher(seqs, gc, TRUE),
            decipher_seqs         = run_decipher(seqs, gc, FALSE)),
            error = function(e) { logm("  %s ERROR: %s", tool, conditionMessage(e))
                                  list(err = TRUE) })
        if (isTRUE(r$err) || !length(r$rows)) {
            record(locus = locus, n_seqs = length(seqs), tool = tool,
                   t_sec = NA, aln_cols = NA, indel_runs = NA,
                   pct_indel_mult3 = NA, codon_phase_purity = NA,
                   status = "error")
            next
        }
        is <- indel_stats(r$rows)
        pp <- phase_purity(r$rows)
        record(locus = locus, n_seqs = length(seqs), tool = tool,
               t_sec = round(r$t, 2),
               aln_cols = nchar(r$rows[1]),
               indel_runs = unname(is["runs"]),
               pct_indel_mult3 = round(unname(is["pct_mult3"]), 4),
               codon_phase_purity = round(pp, 4), status = "ok")
        logm("  %-20s %7.1fs  cols=%d  indel%%3=%.2f  phase=%.3f",
             tool, r$t, nchar(r$rows[1]), is["pct_mult3"], pp)
    }
}
logm("benchmark done -> %s", TSV)
