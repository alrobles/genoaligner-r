#!/usr/bin/env Rscript
# bench_sim.R — axis-A benchmark: tools vs exact simulated ground truth.
#
# For each condition in the grid, generate a dataset with simulate_codon.R,
# run every tool on the unaligned tips, and score the inferred alignment
# against the true one.
#
# Metrics vs truth:
#   sps   — sum-of-pairs score (fraction of true homologous residue pairs
#           recovered; BAliBASE convention)
#   tc    — total-column score (fraction of true columns recovered exactly)
#   fs_recall / fs_precision — non-mult-3 indel events recovered within
#           ±3 columns / inferred non-mult-3 gap runs that match truth
#   n_stops — mean premature stops per row after degap+translate (gc 1)
#   plus: t_sec, width, pct_indel_mult3, codon_phase_purity (structural)
#
# Usage:
#   Rscript bench/bench_sim.R [--grid=small|full] [--tools=a,b,c]

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
tools_keep <- strsplit(opt("tools",
    "msa_codon,msa_codon_refine,msa_codon_protein,msa_align,macse,decipher_translation,decipher_seqs,mafft,prank"),
    ",")[[1]]

# external binaries live in the conda benchtools env
cb <- "~/miniconda3/envs/benchtools/bin"
if (dir.exists(cb)) Sys.setenv(PATH = paste(normalizePath(cb),
                                           Sys.getenv("PATH"), sep = ":"))

BENCH <- normalizePath("bench")
SIM   <- file.path(BENCH, "sim"); dir.create(SIM, showWarnings = FALSE, recursive = TRUE)
MACSE <- file.path(BENCH, "tools", "macse_v2.07.jar")
TSV   <- file.path(BENCH, "out", paste0("sim_results", opt("tag", ""), ".tsv"))
dir.create(dirname(TSV), showWarnings = FALSE, recursive = TRUE)
LOG   <- file.path(BENCH, "out", "sim.log")

logm <- function(...) {
    msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...))
    cat(msg, "\n"); cat(msg, "\n", file = LOG, append = TRUE)
}

## ---- IO helpers -------------------------------------------------------------

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
    for (i in seq_along(v)) writeLines(c(paste0(">", names(v)[i]), v[i]), con)
    close(con)
    invisible(path)
}

# gapped rows -> char matrix, non-ACGT -> "-" (MACSE ! and * markers)
to_mat <- function(rows) {
    m <- do.call(rbind, strsplit(rows, "", fixed = TRUE))
    m[!m %in% c("A", "C", "G", "T")] <- "-"
    m
}

## ---- truth metrics ----------------------------------------------------------

# per-column residue ordinal (NA for gap)
row_ord <- function(m) {
    t(vapply(seq_len(nrow(m)), function(i) {
        o <- cumsum(m[i, ] != "-"); o[m[i, ] == "-"] <- NA; o
    }, numeric(ncol(m))))
}

sps_score <- function(tM, iM) {
    tO <- row_ord(tM); iO <- row_ord(iM)
    hit <- tot <- 0
    for (a in seq_len(nrow(tM) - 1L)) for (b in (a + 1L):nrow(tM)) {
        ct <- which(!is.na(tO[a, ]) & !is.na(tO[b, ]))
        ci <- which(!is.na(iO[a, ]) & !is.na(iO[b, ]))
        kt <- tO[a, ct] * 1e7 + tO[b, ct]
        ki <- iO[a, ci] * 1e7 + iO[b, ci]
        hit <- hit + length(intersect(kt, ki)); tot <- tot + length(kt)
    }
    hit / tot
}

tc_score <- function(tM, iM) {
    sig <- function(m) {
        o <- row_ord(m); o[is.na(o)] <- 0
        apply(o, 2, paste0, collapse = ",")
    }
    st <- sig(tM); si <- sig(iM)
    mean(st %in% si)
}

# gap runs per row (start, end, len)
gap_runs <- function(row) {
    r <- rle(row == "-")
    ends <- cumsum(r$lengths)
    data.frame(start = ends - r$lengths + 1L, end = ends,
               len = r$lengths)[r$values, , drop = FALSE]
}

# frameshift recall/precision for one dataset
# ev: events rows (fs only) with tip, k_anchor, true_col
# iM, iO: inferred matrix + residue ordinals; row order == tip order
fs_metrics <- function(ev, iM, iO) {
    if (!nrow(ev)) return(c(fs_true = 0, fs_rec = 0, fs_recall = NA,
                            fs_runs = NA_real_, fs_spur = NA, fs_prec = NA))
    n_runs <- 0L; n_spur <- 0L; rec <- logical(nrow(ev))
    # map each true event to the inferred column of its anchor residue
    for (k in seq_len(nrow(ev))) {
        i <- ev$tip[k]
        c_inf <- which(iO[i, ] == ev$k_anchor[k])
        if (!length(c_inf)) next
        runs <- gap_runs(iM[i, ]); f <- runs[runs$len %% 3L != 0L, , drop = FALSE]
        rec[k] <- any(f$start <= c_inf + 3 & f$end >= c_inf - 3)
    }
    for (i in unique(ev$tip)) {
        runs <- gap_runs(iM[i, ]); f <- runs[runs$len %% 3L != 0L, , drop = FALSE]
        n_runs <- n_runs + nrow(f)
        if (!nrow(f)) next
        evi <- ev[ev$tip == i, ]
        true_cols <- vapply(evi$k_anchor, function(k) {
            w <- which(iO[i, ] == k); if (length(w)) w else NA_integer_
        }, integer(1))
        true_cols <- true_cols[!is.na(true_cols)]
        spur <- sum(!vapply(seq_len(nrow(f)), function(r)
            any(f$start[r] <= true_cols + 3 & f$end[r] >= true_cols - 3),
            logical(1)))
        n_spur <- n_spur + spur
    }
    c(fs_true = nrow(ev), fs_rec = sum(rec), fs_recall = mean(rec),
      fs_runs = n_runs, fs_spur = n_spur,
      fs_prec = if (n_runs) 1 - n_spur / n_runs else NA)
}

# mean premature stops, standard code
n_stops <- function(rows) {
    cds <- gsub("[^ACGT]", "", rows)
    mean(vapply(cds, function(s) {
        cc <- substring(s, seq(1, nchar(s) - 2, 3), seq(3, nchar(s) - 2, 3))
        sum(cc %in% c("TAA", "TAG", "TGA"))
    }, numeric(1)))
}

## ---- structural metrics (as axis B) -----------------------------------------

indel_stats <- function(rows) {
    lens <- unlist(lapply(rows, function(s) {
        g <- gregexpr("-+", s)[[1]]
        if (g[1] < 0) integer(0) else attr(g, "match.length")
    }))
    if (!length(lens)) return(c(runs = 0, pct_mult3 = NA_real_))
    c(runs = length(lens), pct_mult3 = mean(lens %% 3L == 0L))
}

phase_purity <- function(m) {
    ph <- t(vapply(seq_len(nrow(m)), function(i) {
        isb <- m[i, ] != "-"; (cumsum(isb) - 1L) %% 3L
    }, integer(ncol(m))))
    mean(vapply(seq_len(ncol(m)), function(j) {
        length(unique(ph[m[, j] != "-", j])) <= 1L
    }, logical(1)))
}

## ---- tool runners ------------------------------------------------------------

run_tool <- function(tool, seqs, ifile, ofile) {
    t <- system.time({
        rows <- switch(tool,
            msa_codon = apply(as.matrix(msa_codon(seqs, names = names(seqs),
                              gc = 1L, refine = FALSE)), 1, paste0, collapse = ""),
            msa_codon_refine = apply(as.matrix(msa_codon(seqs, names = names(seqs),
                              gc = 1L, refine = TRUE)), 1, paste0, collapse = ""),
            msa_codon_protein = apply(as.matrix(msa_codon(seqs, names = names(seqs),
                              gc = 1L, backbone = "protein")), 1, paste0, collapse = ""),
            msa_codon_protein_refine = apply(as.matrix(msa_codon(seqs, names = names(seqs),
                              gc = 1L, backbone = "protein", refine = TRUE)),
                              1, paste0, collapse = ""),
            msa_align = apply(as.matrix(msa_align(seqs, names = names(seqs))),
                              1, paste0, collapse = ""),
            decipher_translation = {
                a <- AlignTranslation(DNAStringSet(seqs), verbose = FALSE)
                toupper(as.character(a))
            },
            decipher_seqs = {
                a <- AlignSeqs(DNAStringSet(seqs), verbose = FALSE)
                toupper(as.character(a))
            },
            macse = {
                write_fasta(seqs, ifile)
                system2("java", c("-jar", shQuote(MACSE), "-prog", "alignSequences",
                                  "-seq", shQuote(ifile), "-gc_def", "1",
                                  "-out_NT", shQuote(ofile)),
                        stdout = FALSE, stderr = FALSE)
                if (file.exists(ofile)) read_fasta(ofile) else character(0)
            },
            mafft = {
                write_fasta(seqs, ifile)
                if (!nzchar(Sys.which("mafft"))) character(0) else {
                    system2("mafft", c("--auto", "--quiet", shQuote(ifile)),
                            stdout = ofile, stderr = FALSE)
                    if (file.exists(ofile)) read_fasta(ofile) else character(0)
                }
            },
            prank = {
                write_fasta(seqs, ifile)
                if (!nzchar(Sys.which("prank"))) character(0) else {
                    system2("prank", c(paste0("-d=", ifile),
                                       paste0("-o=", sub("\\.fa$", "", ofile)),
                                       "-codon", "-once"),
                            stdout = FALSE, stderr = FALSE)
                    f <- paste0(sub("\\.fa$", "", ofile), ".best.fas")
                    if (file.exists(f)) read_fasta(f) else character(0)
                }
            },
            character(0))
    })
    list(rows = rows, t = t[["elapsed"]])
}

## ---- driver ------------------------------------------------------------------

grid <- if (opt("grid", "small") == "small") {
    expand.grid(n = c(20, 50), div = c(0.1, 0.3), indel = 0.02,
                pfs = c(0, 0.2), rep = 1:5)
} else {
    expand.grid(n = c(20, 50, 100), div = c(0.1, 0.3), indel = c(0.01, 0.02),
                pfs = c(0, 0.2), rep = 1:10)
}

tools <- intersect(c("msa_codon", "msa_codon_refine", "msa_codon_protein",
                     "msa_codon_protein_refine", "msa_align", "macse",
                     "decipher_translation", "decipher_seqs", "mafft", "prank"),
                   tools_keep)

record <- function(...) {
    row <- data.frame(..., stringsAsFactors = FALSE)
    write.table(row, TSV, sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = !file.exists(TSV), append = file.exists(TSV))
}

logm("sim benchmark start: %d datasets, tools=%s", nrow(grid),
     paste(tools, collapse = ","))

for (g in seq_len(nrow(grid))) {
    p <- grid[g, ]
    seed <- 1000 * g + p$rep
    tag <- sprintf("sim_n%d_L300_d%g_i%g_p%g_r%d", p$n, p$div, p$indel, p$pfs, p$rep)
    infa <- file.path(SIM, paste0(tag, ".in.fa"))
    if (!file.exists(infa))
        system2("Rscript", c(file.path(BENCH, "simulate_codon.R"),
            sprintf("--out=%s", SIM), sprintf("--n=%d", p$n),
            sprintf("--div=%g", p$div), sprintf("--indel=%g", p$indel),
            sprintf("--pfs=%g", p$pfs), sprintf("--rep=%d", p$rep),
            sprintf("--seed=%d", seed)))
    tM <- to_mat(read_fasta(file.path(SIM, paste0(tag, ".true.aln"))))
    ev <- read.delim(file.path(SIM, paste0(tag, ".events.tsv")))
    ev <- ev[ev$fs == TRUE | ev$fs == "TRUE", ]
    seqs <- read_fasta(infa)
    logm("%s  n=%d div=%g pfs=%g rep=%d  (W_true=%d, fs_events=%d)",
         tag, p$n, p$div, p$pfs, p$rep, ncol(tM), nrow(ev))

    for (tool in tools) {
        ifile <- file.path(SIM, paste0(tag, ".", tool, ".in"))
        ofile <- file.path(SIM, paste0(tag, ".", tool, ".out"))
        r <- tryCatch(run_tool(tool, seqs, ifile, ofile),
                      error = function(e) { logm("  %s ERROR: %s", tool,
                          conditionMessage(e)); list(rows = character(0), t = NA) })
        if (!length(r$rows)) {
            record(tag = tag, n = p$n, div = p$div, pfs = p$pfs, rep = p$rep,
                   tool = tool, t_sec = NA, width = NA, sps = NA, tc = NA,
                   fs_recall = NA, fs_prec = NA, fs_true = nrow(ev),
                   n_stops = NA, pct_indel_mult3 = NA, phase = NA, status = "error")
            next
        }
        rr <- if (all(names(seqs) %in% names(r$rows)))
            r$rows[names(seqs)] else r$rows       # align rows to truth order
        iM <- to_mat(rr)
        iO <- row_ord(iM)
        sps <- sps_score(tM, iM); tc <- tc_score(tM, iM); nst <- n_stops(r$rows)
        fsm <- fs_metrics(ev, iM, iO)
        is <- indel_stats(r$rows)
        record(tag = tag, n = p$n, div = p$div, pfs = p$pfs, rep = p$rep,
               tool = tool, t_sec = round(r$t, 2), width = ncol(iM),
               sps = round(sps, 4), tc = round(tc, 4),
               fs_recall = round(unname(fsm["fs_recall"]), 4),
               fs_prec = round(unname(fsm["fs_prec"]), 4),
               fs_true = nrow(ev), n_stops = round(nst, 2),
               pct_indel_mult3 = round(unname(is["pct_mult3"]), 4),
               phase = round(phase_purity(iM), 4), status = "ok")
        logm("  %-20s %7.1fs  sps=%.3f tc=%.3f fsR=%.2f fsP=%.2f stops=%.1f",
             tool, r$t, sps, tc,
             unname(fsm["fs_recall"]), unname(fsm["fs_prec"]), nst)
    }
}
logm("sim benchmark done -> %s", TSV)
