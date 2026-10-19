#!/usr/bin/env Rscript
# make_figures.R — paper figures from bench/out/*.tsv
# Outputs PDFs into paper/figures/. Run from repo root.

suppressPackageStartupMessages({ library(ggplot2); library(gridExtra) })
dir.create("paper/figures", showWarnings = FALSE, recursive = TRUE)

## ---- load axis A (exclude reverted-engine _fix file) -------------------------
fs <- c("bench/out/sim_results.tsv", "bench/out/sim_resultsauto.tsv",
        "bench/out/sim_resultsbbp.tsv")
dl <- lapply(fs[file.exists(fs)], read.delim)
cn <- Reduce(intersect, lapply(dl, names))
sim <- do.call(rbind, lapply(dl, `[`, cn))
sim <- sim[sim$status == "ok" & !is.na(sim$sps), ]

order <- c("msa_codon", "msa_codon_auto", "msa_codon_protein",
           "msa_codon_protein_auto", "msa_align", "macse",
           "decipher_translation", "decipher_seqs", "mafft", "prank")
lbl <- c("msa_codon", "msa_codon (auto)", "msa_codon protein",
         "msa_codon protein (auto)", "msa_align", "MACSE",
         "DECIPHER-tr", "DECIPHER-seqs", "MAFFT", "PRANK")
sim <- sim[sim$tool %in% order, ]
sim$tool <- factor(lbl[match(sim$tool, order)], levels = lbl)
sim$grp <- paste0("d = ", sim$div)
sim$pfs_l <- paste0("frameshift rate pfs = ", sim$pfs)

palette <- c("#999999", "#E69F00", "#0072B2", "#56B4E9",
             "#009E73", "#D55E00", "#CC79A7", "#F0E442",
             "#000000", "#8B8B8B")
names(palette) <- lbl

## ---- Fig 1: SPS vs divergence, faceted by pfs, n = 50 -----------------------
s1 <- sim[sim$n == 50, ]
ag <- aggregate(sps ~ tool + div + pfs_l, s1, function(x)
    c(m = mean(x), se = sd(x) / sqrt(length(x))))
ag <- do.call(data.frame, ag)
names(ag)[4:5] <- c("sps", "se")

fig1 <- ggplot(ag, aes(div, sps, color = tool, group = tool)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = sps - se, ymax = sps + se), width = 0.015) +
    facet_wrap(~ pfs_l) +
    scale_color_manual(values = palette) +
    scale_y_continuous(limits = c(0, 1.02)) +
    labs(x = "tree divergence d (root-to-tip subs/site)",
         y = "sum-of-pairs score (SPS)", color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8),
          legend.key.width = unit(1.4, "lines"))
ggsave("paper/figures/fig1_sps.pdf", fig1, width = 7, height = 4.5)

## ---- Fig 2: axis B — runtime (log) + phase purity ---------------------------
b <- read.delim("bench/out/bench_results.tsv")
b <- b[b$status == "ok" & !is.na(b$t_sec), ]
bl <- c("msa_codon", "msa_align", "macse", "decipher_translation",
        "decipher_seqs", "mafft")
bll <- c("msa_codon", "msa_align", "MACSE", "DECIPHER-tr",
         "DECIPHER-seqs", "MAFFT")
b <- b[b$tool %in% bl, ]
b$tool <- factor(bll[match(b$tool, bl)], levels = bll)
pal2 <- setNames(c("#0072B2", "#009E73", "#D55E00", "#CC79A7",
                   "#F0E442", "#000000"), bll)

p_time <- ggplot(b, aes(tool, t_sec, fill = tool)) +
    geom_boxplot(outlier.size = 0.7, linewidth = 0.3) +
    scale_y_log10(labels = scales::label_number()) +
    scale_fill_manual(values = pal2, guide = "none") +
    labs(x = NULL, y = "time per locus (s, log)") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_phase <- ggplot(b, aes(tool, codon_phase_purity, fill = tool)) +
    geom_boxplot(linewidth = 0.3) +
    scale_fill_manual(values = pal2, guide = "none") +
    labs(x = NULL, y = "codon phase purity") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig2 <- grid.arrange(p_time, p_phase, ncol = 2)
ggsave("paper/figures/fig2_real.pdf", fig2, width = 7.5, height = 3.4)

## ---- Fig 3: full-scale runtime vs n -----------------------------------------
f <- read.delim("bench/out/bench_results_fullscale.tsv")
f <- f[f$tool == "msa_codon" & f$status == "ok", ]
fig3 <- ggplot(f, aes(n_seqs, t_sec)) +
    geom_point(color = "#0072B2", size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                color = "#D55E00", linetype = 2, linewidth = 0.7) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "sequences per locus (log)", y = "msa_codon runtime (s, log)") +
    theme_bw(base_size = 11)
ggsave("paper/figures/fig3_scale.pdf", fig3, width = 4.5, height = 3.2)

cat("figures -> paper/figures/\n")
