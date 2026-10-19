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
#       --indel=0.02 --pfs=0.2 --pcomp=0.5 --pterm=0.3 --rep=1 --seed=1234
#
# --pcomp: probability that a frameshift indel is immediately followed on the
#   same branch by an opposite-sign compensating fs indel 3-15 nt downstream,
#   its length drawn so the net length change is 0 mod 3 (del1+ins1 class).
#   These are invisible to a len%%3 flag -- the reason rowtruth exists.
# --pterm: fraction of tips whose 3' end is truncated by 1-2 nt (a terminal
#   partial codon, NOT an internal frameshift -- the ragged-CDS confound).
# Emits in addition:
#   <tag>.rowtruth.tsv per-tip truth: len, len%%3, retained fs events,
#                        compensated-pair events, terminal trim, internal_fs
#   <tag>.tree.nwk     the true guide tree (for guide-tree diagnostics)

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
p_comp <- as.numeric(kv("pcomp", 0))          # P(fs indel gets compensating partner)
p_term <- as.numeric(kv("pterm", 0))          # fraction of tips 3'-trimmed 1-2 nt
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
#   (accumulates down the tree: a tip's ev holds every ancestral event)
# pc: pair-counter environment (pc$n) so compensated partners share an id
# across the whole simulation.

FS_LEN <- c(1L, 2L, 4L, 5L, 7L, 8L, 10L)   # internal-frameshift lengths
FS_W   <- 2^-(1:7)

# Compensating partner length: same mod-3 residue as len1 so the net
# length change is 0 mod 3 (del1+ins1, ins4+del7, ...) -- invisible to a
# len%%3 flag but real internal frameshifting.
comp_len <- function(len1)
    sample(FS_LEN[FS_LEN %% 3L == len1 %% 3L], 1L)

evolve_branch <- function(s, cid, bl, indel_rate, p_fs, p_comp, ev, pc) {
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
        len <- if (fs) sample(FS_LEN, 1, prob = FS_W)
               else 3L * (rgeom(1, 0.4) + 1L)
        pid <- 0L
        if (del) {
            end <- min(pos + len - 1L, L)
            anchor <- if (pos > 1) cid[pos - 1] else if (end < L) cid[end + 1] else NA_real_
            keep <- setdiff(seq_len(L), pos:end)
            s <- s[keep]; cid <- cid[keep]
            ev <- rbind(ev, data.frame(anchor = anchor, type = "del",
                                       len = end - pos + 1L, fs = fs,
                                       pair_id = 0L))
        } else {
            left <- cid[pos]; right <- if (pos < L) cid[pos + 1] else left + 1
            newid <- left + (right - left) * seq_len(len) / (len + 1)
            s <- append(s, sample(DNA, len, TRUE), after = pos)
            cid <- append(cid, newid, after = pos)
            ev <- rbind(ev, data.frame(anchor = left, type = "ins",
                                       len = len, fs = fs, pair_id = 0L))
        }
        # Compensating partner: opposite sign, 3-15 nt downstream on the
        # same branch, length drawn to cancel the frame offset. Applied to
        # the already-mutated sequence; cid bookkeeping is identical.
        if (fs && runif(1) < p_comp) {
            L <- length(s); if (L >= 30) {
                pid <- pc$n <- pc$n + 1L
                ev$pair_id[nrow(ev)] <- pid
                len2 <- comp_len(len)
                off  <- sample(3:15, 1)
                q    <- min(pos + off, L)          # partner site, post-event coords
                if (del) {                          # del, then ins
                    left <- cid[q]; right <- if (q < L) cid[q + 1] else left + 1
                    newid <- left + (right - left) * seq_len(len2) / (len2 + 1)
                    s <- append(s, sample(DNA, len2, TRUE), after = q)
                    cid <- append(cid, newid, after = q)
                    ev <- rbind(ev, data.frame(anchor = left, type = "ins",
                                               len = len2, fs = TRUE,
                                               pair_id = pid))
                } else {                            # ins, then del
                    end2 <- min(q + len2 - 1L, L)
                    anchor <- if (q > 1) cid[q - 1] else if (end2 < L) cid[end2 + 1] else NA_real_
                    keep <- setdiff(seq_len(L), q:end2)
                    s <- s[keep]; cid <- cid[keep]
                    ev <- rbind(ev, data.frame(anchor = anchor, type = "del",
                                               len = end2 - q + 1L, fs = TRUE,
                                               pair_id = pid))
                }
            }
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
                       len = integer(0), fs = logical(0),
                       pair_id = integer(0))
pc <- new.env(parent = emptyenv()); pc$n <- 0L

tipenv <- new.env(parent = emptyenv())
rec <- function(node, s, cid, ev) {
    for (e in which(tr$edge[, 1] == node)) {
        nd <- tr$edge[e, 2]
        r <- evolve_branch(s, cid, tr$edge.length[e], indel_rate, p_fs,
                           p_comp, ev, pc)
        if (nd <= n_taxa) tipenv[[as.character(nd)]] <- r
        else rec(nd, r$s, r$cid, r$ev)
    }
}
rec(n_taxa + 1L, root_s, root_cid, empty_ev)
tips <- lapply(as.character(seq_len(n_taxa)), function(k) tipenv[[k]])

# Terminal partials: a fraction of tips lose 1-2 nt off the 3' end --
# the ragged-CDS confound (len%%3 != 0 WITHOUT an internal frameshift).
# Not recorded as an event: it is truncation, not an indel on the tree.
term_trim <- integer(n_taxa)
if (p_term > 0) {
    for (i in which(runif(n_taxa) < p_term)) {
        tr_n <- sample(1:2, 1)
        ti <- tips[[i]]
        if (length(ti$s) > tr_n + 30L) {
            term_trim[i] <- tr_n
            tips[[i]]$s   <- head(ti$s,   -tr_n)
            tips[[i]]$cid <- head(ti$cid, -tr_n)
        }
    }
}

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
    pair_id = integer(0), k_anchor = integer(0), true_col = integer(0))

# Per-tip truth for flag diagnostics: internal_fs = the tip's lineage
# carried >=1 retained frameshift event (all sim fs events are internal);
# compensated = >=1 event belonging to a compensated pair; term_trim is
# the 3'-truncation size (a terminal partial, never an internal fs).
rt <- data.frame(
    tip = seq_len(n_taxa), label = tr$tip.label,
    len = nchar(seqs), len_mod3 = nchar(seqs) %% 3L,
    n_fs = 0L, n_comp = 0L, term_trim,
    internal_fs = FALSE, compensated = FALSE)
if (nrow(events)) {
    fsev <- events[events$fs, , drop = FALSE]
    if (nrow(fsev)) {
        tab <- table(fsev$tip)
        rt$n_fs[as.integer(names(tab))] <- as.integer(tab)
        cev <- fsev[fsev$pair_id > 0L, , drop = FALSE]
        if (nrow(cev)) {
            ctab <- table(cev$tip)
            rt$n_comp[as.integer(names(ctab))] <- as.integer(ctab)
        }
    }
}
rt$internal_fs <- rt$n_fs > 0L
rt$compensated <- rt$n_comp > 0L

tag <- sprintf("sim_n%d_L%d_d%g_i%g_p%g%s%s_r%d", n_taxa, n_codons,
               divergence, indel_rate, p_fs,
               if (p_comp > 0) sprintf("_c%g", p_comp) else "",
               if (p_term > 0) sprintf("_t%g", p_term) else "",
               rep_i)

write_fasta <- function(x, path, mat = FALSE) {
    nm <- if (mat) rownames(x) else names(x)
    if (mat) x <- apply(x, 1, paste0, collapse = "")
    writeLines(paste(rbind(paste0(">", nm), x), collapse = "\n"), path)
}
write_fasta(seqs, file.path(outdir, paste0(tag, ".in.fa")))
write_fasta(true_aln, file.path(outdir, paste0(tag, ".true.aln")), mat = TRUE)
write.table(events[, c("tip", "k_anchor", "true_col", "type", "len", "fs",
                       "pair_id")],
            file.path(outdir, paste0(tag, ".events.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rt, file.path(outdir, paste0(tag, ".rowtruth.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.tree(tr, file.path(outdir, paste0(tag, ".tree.nwk")))
writeLines(c(
    sprintf("n_taxa\t%d", n_taxa), sprintf("n_codons\t%d", n_codons),
    sprintf("divergence\t%g", divergence), sprintf("indel_rate\t%g", indel_rate),
    sprintf("p_fs\t%g", p_fs), sprintf("p_comp\t%g", p_comp),
    sprintf("p_term\t%g", p_term),
    sprintf("rep\t%d", rep_i), sprintf("seed\t%d", seed),
    sprintf("true_width\t%d", W),
    sprintf("n_events\t%d", nrow(events)),
    sprintf("n_fs_events\t%d", sum(events$fs)),
    sprintf("n_comp_events\t%d", sum(events$pair_id > 0L)),
    sprintf("n_tips_internal_fs\t%d", sum(rt$internal_fs)),
    sprintf("n_tips_term_trim\t%d", sum(rt$term_trim > 0L))
), file.path(outdir, paste0(tag, ".meta.txt")))

cat(sprintf("[%s] %s  n=%d W=%d events=%d fs=%d comp=%d ifs=%d trim=%d\n",
            format(Sys.time(), "%H:%M:%S"), tag, n_taxa, W, nrow(events),
            sum(events$fs), sum(events$pair_id > 0L),
            sum(rt$internal_fs), sum(rt$term_trim > 0L)))
