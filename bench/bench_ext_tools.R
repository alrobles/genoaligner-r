#!/usr/bin/env Rscript
# Extend the n=50 head-to-head benchmark with MAFFT and PRANK-codon.
# Uses the SAME seeded subsets (bench/out/subsets/*.in.fa) that
# bench_codon_msa.R wrote, so results are directly comparable.
# Appends rows to bench/out/bench_results.tsv with the same schema.

BENCH <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
OUT   <- file.path(BENCH, "out")
SUB   <- file.path(OUT, "subsets")
TSV   <- file.path(OUT, "bench_results.tsv")

logm <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
                          sprintf(...), "\n")

read_fasta <- function(path) {
    l <- readLines(path); h <- grep("^>", l)
    if (!length(h)) return(character(0))
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

record <- function(...) {
    row <- data.frame(..., stringsAsFactors = FALSE)
    write.table(row, TSV, sep = "\t", quote = FALSE, row.names = FALSE,
                col.names = !file.exists(TSV), append = file.exists(TSV))
}

run_tool <- function(tool, ifile, ofile) {
    t <- system.time({
        rows <- switch(tool,
            mafft = {
                system2("mafft", c("--auto", "--quiet", shQuote(ifile)),
                        stdout = ofile, stderr = FALSE)
                if (file.exists(ofile)) read_fasta(ofile) else character(0)
            },
            prank = {
                system2("prank", c(paste0("-d=", ifile),
                                   paste0("-o=", sub("\\.fa$", "", ofile)),
                                   "-codon", "-once"),
                        stdout = FALSE, stderr = FALSE)
                f <- paste0(sub("\\.fa$", "", ofile), ".best.fas")
                if (file.exists(f)) read_fasta(f) else character(0)
            },
            character(0))
    })
    list(rows = rows, t = t[["elapsed"]])
}

inputs <- sort(list.files(SUB, pattern = "\\.in\\.fa$"))
inputs <- inputs[!file.exists(file.path(SUB, sub("\\.in\\.fa$", ".done.flag", inputs)))]
done <- if (file.exists(TSV)) unique(read.delim(TSV)$locus[
    read.delim(TSV)$tool %in% c("mafft", "prank")]) else character(0)
tools <- c("mafft", "prank")

for (f in inputs) {
    locus <- sub("\\.in\\.fa$", "", f)
    ifile <- file.path(SUB, f)
    n <- length(read_fasta(ifile))
    for (tool in tools) {
        if (file.exists(TSV) &&
            any(dd <- (read.delim(TSV)$locus == locus &
                       read.delim(TSV)$tool == tool))) next
        ofile <- file.path(SUB, paste0(locus, ".", tool, ".fa"))
        r <- tryCatch(run_tool(tool, ifile, ofile),
                      error = function(e) list(rows = character(0), t = NA))
        if (!length(r$rows)) {
            record(locus = locus, n_seqs = n, tool = tool, t_sec = NA,
                   aln_cols = NA, indel_runs = NA, pct_indel_mult3 = NA,
                   codon_phase_purity = NA, status = "error")
            logm("%s %s: no output", locus, tool); next
        }
        is <- indel_stats(r$rows); pp <- phase_purity(r$rows)
        record(locus = locus, n_seqs = n, tool = tool,
               t_sec = round(r$t, 2), aln_cols = nchar(r$rows[1]),
               indel_runs = unname(is["runs"]),
               pct_indel_mult3 = round(unname(is["pct_mult3"]), 4),
               codon_phase_purity = round(pp, 4), status = "ok")
        logm("%s %-6s %7.1fs cols=%d indel%%3=%.2f phase=%.3f",
             locus, tool, r$t, nchar(r$rows[1]), is["pct_mult3"], pp)
    }
}
logm("done -> %s", TSV)
