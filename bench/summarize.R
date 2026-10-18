#!/usr/bin/env Rscript
# Summarize bench/out/bench_results*.tsv into headline tables.
# Usage: Rscript bench/summarize.R [tag]   (default: all bench_results*.tsv)

args <- commandArgs(TRUE)
tag <- if (length(args)) args[1] else ""
files <- list.files("bench/out", pattern = "bench_results.*\\.tsv$",
                    full.names = TRUE)
if (nzchar(tag)) files <- files[grepl(tag, files)]

d <- do.call(rbind, lapply(files, function(f)
    tryCatch(read.delim(f, stringsAsFactors = FALSE),
             error = function(e) NULL)))
d <- d[d$status == "ok", ]
if (!nrow(d)) stop("no ok rows yet")

cat("== rows:", nrow(d), "| loci:", length(unique(d$locus)), "\n\n")

## per-tool summary across loci
agg <- aggregate(cbind(t_sec, indel_runs) ~ tool, d, sum, na.rm = TRUE)
agg2 <- aggregate(cbind(pct_indel_mult3, codon_phase_purity) ~ tool, d,
                  mean, na.rm = TRUE)
s <- merge(agg, agg2, by = "tool")
s$median_t <- aggregate(t_sec ~ tool, d, median)$t_sec
s$n_loci <- aggregate(locus ~ tool, d, function(x) length(unique(x)))$locus
s <- s[order(s$tool), ]
cat("== per-tool (totals over loci; metrics averaged) ==\n")
print(s, row.names = FALSE, digits = 3)

## MACSE vs msa_codon speedup per locus (same n)
if (all(c("macse", "msa_codon") %in% d$tool)) {
    m <- merge(d[d$tool == "macse", c("locus", "n_seqs", "t_sec")],
               d[d$tool == "msa_codon", c("locus", "t_sec")],
               by = "locus", suffixes = c("_macse", "_codon"))
    m$speedup <- round(m$t_sec_macse / m$t_sec_codon, 1)
    cat("\n== MACSE vs msa_codon per-locus ==\n")
    print(m[order(-m$speedup), ], row.names = FALSE)
    cat("\nmedian speedup:", median(m$speedup), "x\n")
}

## per-locus table (wide)
w <- reshape(d[, c("locus", "tool", "t_sec", "codon_phase_purity")],
             timevar = "tool", idvar = "locus", direction = "wide")
cat("\n== per-locus ==\n")
print(w[order(w$locus), ], row.names = FALSE)
