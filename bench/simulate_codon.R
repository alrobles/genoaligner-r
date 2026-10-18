#!/usr/bin/env Rscript
# simulate_codon.R — codon-MSA benchmark simulator with exact ground truth.
#
# Evolves a coding sequence along a random Yule tree. Every site carries a
# stable column id; insertions get fractional ids so that sorting recovers
# the true column order. Emits:
#   <tag>.in.fa      unaligned tip sequences (tool input)
#   <tag>.true.aln   true multiple alignment (gapped fasta)
#   <tag>.events.tsv indel events per tip (anchor residue ordinal, true col)
#   <tag>.meta.txt   parameters
#
# Usage:
#   Rscript simulate_codon.R --out=DIR --n=50 --codons=300 --div=0.1 \
#       --indel=0.02 --pfs=0.2 --rep=1 --seed=1234

suppressMessages(library(ape))

DNA <- c("A", "C", "G", "T")
STOPS <- c("TAA", "TAG", "TGA")

args <- commandArgs(trailingOnly = TRUE)
kv <- function(k, d = NA) {
    m <- grep(paste0("^--", k, "="), args, value = TRUE)
    if (length(m)) sub(paste0("^--", k, "="), "", m[1]) else d
}
outdir <- kv("out", "bench/sim"); dir.create(outdir, TRUE, TRUE)
n_taxa <- as.integer(kv("n", 50))
n_codons <- as.integer(kv("codons", 300))
divergence <- as.numeric(kv("div", 0.1))      # mean root-to-tip subs/site
indel_rate <- as.numeric(kv("indel", 0.02))   # events per site per unit bl
p_fs <- as.numeric(kv("pfs", 0.2))            # fraction of indels breaking frame
rep_i <- as.integer(kv("rep", 1))
seed <- as.integer(kv("seed", 1))

## ---- root sequence: random codons, no stops in frame ----
mk_root <- function(n_codons) {
    s <- sample(DNA, 3L * n_codons, replace = TRUE)
    cod <- matrix(s, nrow = 3)
    for (j in which(apply(cod, 2, function(x)
        paste0(x, collapse = "") %in% STOPS))) {
        repeat {
            c3 <- sample(DNA, 3, TRUE)
            if (!paste0(c3, collapse = "") %in% STOPS) {
                s[(3L * j - 2L):(3L * j)] <- c3; break
            }
        }
    }
    s
}

## ---- evolve along one branch ----
# s: nucleotide char vector; cid: column id per site; ev: events data.frame
evolve_branch <- function(s, cid, bl, indel_rate, p_fs, ev) {
    L <- length(s)
    # point substitutions (JC69)
    p_sub <- 1 - exp(-bl)
    idx <- which(runif(L) < p_sub)
    for (k in idx) s[k] <- sample(setdiff(DNA, s[k]), 1L)
    # indels
    n_ev <- rpois(1, indel_rate * bl * L)
    for (e in seq_len(n_ev)) {
        L <- length(s); if (L < 30) break
        pos <- sample.int(L, 1)
        del <- runif(1) < 0.5
        fs <- runif(1) < p_fs
        len <- if (fs) {
            sample(c(1, 2, 4, 5, 7, 8, 10), 1, prob = 2^-(1:7))
        } else 3L * (rgeom(1, 0.4) + 1L)
        if (del) {
            end <- min(pos + len - 1L, L)
            anchor <- if (pos > 1) cid[pos - 1] else if (end < L) cid[end + 1] else NA_real_
            keep <- setdiff(seq_len(L), pos:end)
            s <- s[keep]; cid <- cid[keep]
            ev <- rbind(ev, data.frame(anchor = anchor, type = "del",
                                       len = end - pos + 1L, fs = fs))
        } else {
            left <- cid[pos]; right <- if (pos < L) cid[pos + 1] else left + 1
            newid <- left + (right - left) * seq_len(len) / (len + 1)
            s <- append(s, sample(DNA, len, TRUE), after = pos)
            cid <- append(cid, newid, after = pos)
            ev <- rbind(ev, data.frame(anchor = left, type = "ins",
                                       len = len, fs = fs))
        }
    }
    list(s = s, cid = cid, ev = ev)
}

## ---- simulate ----
set.seed(seed)
tr <- rtree(n_taxa)
dep <- node.depth.edgelength(tr)
tr$edge.length <- tr$edge.length * (divergence / mean(dep[seq_len(n_taxa)]))

root_s <- mk_root(n_codons)
root_cid <- seq_len(3L * n_codons)
empty_ev <- data.frame(anchor = numeric(0), type = character(0),
                       len = integer(0), fs = logical(0))

tipenv <- new.env(parent = emptyenv())
rec <- function(node, s, cid, ev) {
    for (e in which(tr$edge[, 1] == node)) {
        nd <- tr$edge[e, 2]
        r <- evolve_branch(s, cid, tr$edge.length[e], indel_rate, p_fs, ev)
        if (nd <= n_taxa) tipenv[[as.character(nd)]] <- r
        else rec(nd, r$s, r$cid, r$ev)
    }
}
rec(n_taxa + 1L, root_s, root_cid, empty_ev)
tips <- lapply(as.character(seq_len(n_taxa)), function(k) tipenv[[k]])

## ---- emit ----
all_colid <- sort(unique(unlist(lapply(tips, `[[`, "cid"))))
W <- length(all_colid)

true_aln <- matrix("-", n_taxa, W)
seqs <- character(n_taxa)
ev_rows <- list()
for (i in seq_len(n_taxa)) {
    ti <- tips[[i]]
    seqs[i] <- paste0(ti$s, collapse = "")
    true_aln[i, match(ti$cid, all_colid)] <- ti$s
    ev <- ti$ev
    keep <- !is.na(ev$anchor) & ev$anchor %in% ti$cid
    ev <- ev[keep, , drop = FALSE]
    if (nrow(ev)) {
        ev$tip <- i
        ev$k_anchor <- match(ev$anchor, ti$cid)   # residue ordinal in ungapped seq
        ev$true_col <- match(ev$anchor, all_colid)
        ev_rows[[i]] <- ev
    }
}
names(seqs) <- tr$tip.label
rownames(true_aln) <- tr$tip.label
events <- do.call(rbind, ev_rows)
if (is.null(events)) events <- data.frame(tip = integer(0), anchor = numeric(0),
    type = character(0), len = integer(0), fs = logical(0),
    k_anchor = integer(0), true_col = integer(0))

tag <- sprintf("sim_n%d_L%d_d%g_i%g_p%g_r%d", n_taxa, n_codons,
               divergence, indel_rate, p_fs, rep_i)

write_fasta <- function(x, path, mat = FALSE) {
    nm <- if (mat) rownames(x) else names(x)
    if (mat) x <- apply(x, 1, paste0, collapse = "")
    writeLines(paste(rbind(paste0(">", nm), x), collapse = "\n"), path)
}
write_fasta(seqs, file.path(outdir, paste0(tag, ".in.fa")))
write_fasta(true_aln, file.path(outdir, paste0(tag, ".true.aln")), mat = TRUE)
write.table(events[, c("tip", "k_anchor", "true_col", "type", "len", "fs")],
            file.path(outdir, paste0(tag, ".events.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
    sprintf("n_taxa\t%d", n_taxa), sprintf("n_codons\t%d", n_codons),
    sprintf("divergence\t%g", divergence), sprintf("indel_rate\t%g", indel_rate),
    sprintf("p_fs\t%g", p_fs), sprintf("rep\t%d", rep_i), sprintf("seed\t%d", seed),
    sprintf("true_width\t%d", W),
    sprintf("n_events\t%d", nrow(events)),
    sprintf("n_fs_events\t%d", sum(events$fs))
), file.path(outdir, paste0(tag, ".meta.txt")))

cat(sprintf("[%s] %s  n=%d W=%d events=%d fs=%d\n", format(Sys.time(), "%H:%M:%S"),
            tag, n_taxa, W, nrow(events), sum(events$fs)))
