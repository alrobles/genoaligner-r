# bench — codon-MSA benchmark (Phase 5)

Benchmark scaffolding for the genoaligner app note: `msa_codon` /
`msa_align` vs **MACSE v2.07** vs **DECIPHER** (`AlignSeqs`,
`AlignTranslation`) on the 31 phylogenyAI "classic loci" (real mammalian
coding-gene sets, Upham et al. backbone).

This directory is **excluded from the package tarball** (`^bench$` in
`.Rbuildignore`). External tools are invoked, never shipped: MACSE runs
from `tools/macse_v2.07.jar` (sha256
`96873a14…3963`, bioconda-verified), DECIPHER from Bioconductor.

## Layout

- `data/genes_final/*.fasta` — 31 loci pulled from
  `/beegfs/a474r867/phylogenyAI/data/genes_final/` (KU HPC; source of the
  Upham classic-locus sequence sets).
- `tools/macse_v2.07.jar` — MACSE v2.07 (Ranwez et al. 2020, MBE
  35:2582–2584).
- `bench_codon_msa.R` — driver; samples `n` seqs/locus (same subset for
  every tool, `set.seed(42)`), times each aligner, computes
  frame-conservation metrics, appends to `out/bench_results.tsv`.
- `out/` — results TSV + log + tool input/output fastas.

## Run

```bash
Rscript bench/bench_codon_msa.R --n=50          # 50 seqs/locus (~hours, MACSE-bound)
Rscript bench/bench_codon_msa.R --n=0           # full loci (genoaligner only — see below)
Rscript bench/bench_codon_msa.R --loci='^CYTB$' # single locus
```

Requires: `genoaligner` (installed), `Biostrings` + `DECIPHER`
(Bioconductor), `java` on PATH.

## Metrics

- `t_sec` — wall time, single run, single thread.
- `indel_runs`, `pct_indel_mult3` — per-sequence gap runs; fraction whose
  length is a multiple of 3 (whole-codon indels preserve the frame).
- `codon_phase_purity` — fraction of columns in which every non-gap base
  occupies the same codon phase (pos-1/2/3) within its own reading frame.
  Codon-aware aligners score ~1.0; nucleotide aligners mix phases at every
  non-multiple-of-3 indel.

## Notes

- Mitochondrial loci (COI, CYTB, ND1, ND2) run under NCBI genetic code 2
  (`gc_def 2` in MACSE, `gc = 2` in `msa_codon`); the rest under code 1.
- MACSE output marks frameshifts with `!` — treated as occupying a codon
  position (non-gap) for the metrics.
- Head-to-head uses identical inputs per locus (same seeded subsample).
- Full-locus runs (`--n=0`) are meaningful for `msa_*` (seconds–minutes);
  MACSE does not finish at 2–3.5k seqs on a workstation timescale — that
  gap is itself a headline result, so cap MACSE runs when comparing.
