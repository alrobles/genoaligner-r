#!/usr/bin/env Rscript
# diag_phase0.R — Phase-0 diagnostics for the codon-MSA roadmap.
#
# Three questions, all measurable only because of the new machinery:
#   R2: does the TRUE guide tree improve true-SPS over the NJ guide?
#       (if not, the residual gap is scoring/gap-model, not ordering)
#   R3: does the model's own SP objective rank variants like true-SPS?
#   R1: flag audit — len%%3 vs rowtruth (internal_fs / term_trim /
#       compensated), i.e. how much does the current flag miss/false-fire?
#
# Usage: Rscript bench/diag_phase0.R [--reps=3] [--n=30]

suppressMessages({library(genoaligner); library(ape)})

args <- commandArgs(trailingOnly = TRUE)
kv <- function(k, d) { m <- grep(paste0("^--", k, "="), args, value = TRUE)
    if (length(m)) as.numeric(sub(".*=", "", m[1])) else d }
NREP <- as.integer(kv("reps", 3)); N <- as.integer(kv("n", 30))

SIM <- "bench/sim_diag"; dir.create(SIM, TRUE, TRUE)

## ---- helpers (kept in sync with bench_sim.R) -------------------------------
read_fasta <- function(path) {
    ln <- readLines(path, warn = FALSE)
    hdr <- grepl("^>", ln)
    nm <- sub("^>\\s*", "", ln[hdr])
    seqs <- vapply(split(ln[!hdr], cumsum(hdr)[!hdr]),
                   paste0, character(1), collapse = "")
    names(seqs) <- sub("\\s.*$", "", nm)
    toupper(seqs)
}
to_mat <- function(rows) {
    m <- do.call(rbind, strsplit(rows, "", fixed = TRUE))
    m[!m %in% c("A", "C", "G", "T")] <- "-"
    m
}
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

CFG <- list(
    A = list(pfs = 0,   pcomp = 0,   pterm = 0),
    B = list(pfs = 0.2, pcomp = 0.5, pterm = 0),
    C = list(pfs = 0,   pcomp = 0,   pterm = 0.4),
    D = list(pfs = 0.2, pcomp = 0.5, pterm = 0.4))

rows <- list()
for (d in c(0.1, 0.3)) for (cfg in names(CFG)) for (rep in seq_len(NREP)) {
    p <- CFG[[cfg]]
    cmd <- sprintf(
        "Rscript bench/simulate_codon.R --out=%s --n=%d --codons=200 --div=%g --indel=0.02 --pfs=%g --pcomp=%g --pterm=%g --rep=%d --seed=%d",
        SIM, N, d, p$pfs, p$pcomp, p$pterm, rep, 1000 + rep)
    system(cmd, ignore.stdout = TRUE)
    tag <- sprintf("sim_n%d_L200_d%g_i%g_p%g%s%s_r%d", N, d, 0.02, p$pfs,
                   if (p$pcomp > 0) sprintf("_c%g", p$pcomp) else "",
                   if (p$pterm > 0) sprintf("_t%g", p$pterm) else "", rep)
    base <- file.path(SIM, tag)
    if (!file.exists(paste0(base, ".in.fa"))) { cat("MISSING", base, "\n"); next }
    sq <- read_fasta(paste0(base, ".in.fa"))
    tr <- read.tree(paste0(base, ".tree.nwk"))
    tM <- to_mat(read_fasta(paste0(base, ".true.aln")))
    rt <- read.delim(paste0(base, ".rowtruth.tsv"))
    sq <- sq[tr$tip.label]                       # tree order == tip order

    aln_of <- function(...) apply(as.matrix(msa_codon(sq, ...)), 1,
                                  paste0, collapse = "")
    ## R2: NJ vs true tree, protein backbone
    a_nj <- aln_of(backbone = "protein")
    a_tr <- aln_of(backbone = "protein", guide_tree = tr)
    i_nj <- to_mat(a_nj[names(sq)]); i_tr <- to_mat(a_tr[names(sq)])
    s_nj <- sps_score(tM, i_nj);     s_tr <- sps_score(tM, i_tr)
    ## R3: model SP vs true SPS, same two candidates
    sp_nj <- msa_sp_score(i_nj, "dna"); sp_tr <- msa_sp_score(i_tr, "dna")
    ## R1: flag audit
    flag <- rt$len_mod3 != 0L
    tp <- sum(flag & rt$internal_fs); fp <- sum(flag & !rt$internal_fs)
    fn <- sum(!flag & rt$internal_fs); tn <- sum(!flag & !rt$internal_fs)
    rows[[length(rows) + 1L]] <- data.frame(
        d = d, cfg = cfg, rep = rep, tag = tag,
        sps_nj = s_nj, sps_tr = s_tr, dsps = s_tr - s_nj,
        sp_nj = sp_nj, sp_tr = sp_tr,
        sp_prefers_tr = sp_tr > sp_nj, sps_prefers_tr = s_tr > s_nj,
        n_rows = N, n_ifs = sum(rt$internal_fs),
        n_trim = sum(rt$term_trim > 0), n_comp = sum(rt$compensated),
        flag_prec = tp / max(tp + fp, 1), flag_rec = tp / max(tp + fn, 1),
        flag_missed = fn, flag_false = fp)
    cat(sprintf(paste0("%-3s d=%.1f rep=%d  SPS nj=%.3f tree=%.3f (d=%+.4f)  ",
                       "SP: nj=%.0f tree=%.0f  flag P/R=%.2f/%.2f miss=%d fa=%d\n"),
                cfg, d, rep, s_nj, s_tr, s_tr - s_nj, sp_nj, sp_tr,
                tp / max(tp + fp, 1), tp / max(tp + fn, 1), fn, fp))
}

R <- do.call(rbind, rows)
write.table(R, "bench/out/diag_phase0.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("\n=== summary by cfg x d ===\n")
ag <- aggregate(cbind(sps_nj, sps_tr, dsps, flag_prec, flag_rec,
                      flag_missed, flag_false) ~ cfg + d, R, mean)
print(ag, digits = 3)
cat("\ntree wins:", sum(R$dsps > 1e-4), "/", nrow(R),
    " ties:", sum(abs(R$dsps) <= 1e-4), " loses:", sum(R$dsps < -1e-4), "\n")
cat("SP/SPS rank agreement:", mean(R$sp_prefers_tr == R$sps_prefers_tr), "\n")
