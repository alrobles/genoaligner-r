# Benchmark plan — genoaligner codon-aware MSA

Paper-grade evaluation design. Three axes: (A) simulation with ground truth,
(B) real-data scaling + structure, (C) case studies. Designed so each axis
answers a different reviewer question.

## A. Simulation with ground truth

**Question answered:** is the inferred alignment *correct*, not merely
frame-consistent? `phase = 1.0` on real data is partially by construction
(codon mode forces whole-codon indels); simulation supplies unbiased
accuracy numbers.

**Simulator** (`simulate_codon.R`): custom, dependency-free (base R + ape
for the guide tree). A random Yule tree (`ape::rtree`) is scaled so that
expected nucleotide divergence matches the condition. Along each branch:

- point substitutions: JC69 per site, rate × branch length;
- indels: Poisson(indel_rate × bl × n_codons) events per branch; insertion
  or deletion (50/50); length ~ geometric(mean ≈ 1 codon); with probability
  `p_fs` the length is drawn non-multiple-of-3 (a true frameshift);
- every site carries a stable column id (insertions get fractional ids so
  sorting recovers true column order) → the true MSA is exact, not guessed.

Emitted per dataset: input FASTA, true gapped alignment, indel-event table
(tip, column, type, length, is_frameshift), parameter meta.

**Metrics vs truth** (per dataset, then averaged over replicates):

- **SPS** — sum-of-pairs score: fraction of homologous residue pairs in the
  reference alignment recovered in the inferred one (BAliBASE convention).
- **TC** — total-column score: fraction of reference columns recovered
  exactly.
- **Frameshift recall/precision** — a true non-mult-3 indel event is
  recalled if the inferred row has a non-mult-3 gap run within ±3 columns
  of the true event position; inferred non-mult-3 runs not matching any
  true event are spurious (precision).
- **n_stops** — premature stop codons per row after degapping+translating.
- runtime, alignment width, phase purity, %indel runs of length ≡ 0 (mod 3).

**Condition grid** (kept small so MACSE finishes):

- n_taxa ∈ {20, 50, 100}; n_codons = 300;
- divergence ∈ {0.1, 0.3} expected subs/site;
- indel_rate = 0.02 events/codon/unit-bl (calibrated);
- p_frameshift ∈ {0, 0.2};
- replicates = 5 per condition (seeded).

**Tools:** msa_codon (refine=FALSE and TRUE), msa_align, MACSE v2.07,
DECIPHER (translation + sequence), MAFFT --auto, PRANK -codon.
MACSE may be restricted to n ≤ 50 if wall-clock demands it.

## B. Real data — scaling + structure

**Question answered:** does the method hold frame on real data, and does it
scale to dataset sizes where alternatives fail?

- 31 classic Upham loci (genes_final, phylogenyAI) — the reference set for
  mammalian phylogenomics.
- Tiers: n = 50 (all tools), n = 200 (MACSE-feasible), full locus
  (msa_align/msa_codon only; MACSE/DECIPHER reported as DNF where they
  cannot complete).
- Genetic codes: NCBI 2 for COI, CYTB, ND1, ND2; 1 otherwise.
- Metrics: wall-clock (3 replicates for timing), phase purity, %indel
  runs ≡ 0 (mod 3), n_stops, alignment width.
- MACSE output parsed permissively (`!` frameshift, `*` stop markers kept).
- Subsampling seeded; all command lines recorded in bench.log.

## C. Case studies

**Question answered:** does the alignment matter downstream?

- Pseudogene/frameshift set — replicate MACSE's kiwi-style demo: mixed
  functional + pseudogene sequences, `msa_codon(refine=TRUE)` must recover
  the frameshifts (stage-2 pass).
- Optional anchor result: dN/dS (codeml/HyPhy BUSTED) on msa_codon vs
  msa_align output for the same locus — frame-breaking is expected to
  inflate dN/dS. If confirmed, this is the abstract line.

## Honesty rules

- Never claim `phase = 1` as evidence of quality for msa_codon default —
  it is structural. Accuracy claims come only from axis A.
- "MACSE-class" means: whole-codon indels by default, optional frameshift
  pass in refine mode — not algorithmic equivalence.
- Report MACSE DNF honestly as feasibility result, not quality result.
- All randomness seeded; environment recorded (R, Java, package versions).

## Status / files

- `bench_codon_msa.R` — axis B runner (n=50 running, full-scale done).
- `summarize.R` — TSV → per-locus × tool tables.
- `simulate_codon.R` — axis A simulator.
- `bench_sim.R` — axis A runner + truth metrics.
- `out/` — results (gitignored). `data/`, `tools/` — gitignored inputs.

## Algorithmic branch: protein backbone (2026-10)

Axis A exposed that the codon-token backbone collapses at high divergence
(n=50, d=0.3: SPS ≈ 0.5–0.7 vs msa_align ≈ 0.9) and cannot represent true
frameshifts (fs_recall = 0 by construction). Experiments ruled out the
guide tree (nt k-mer guide → no gain), gap cost, and iter_refine. The
working fix is algorithmic: `msa_codon(backbone = "protein")` — translate
in frame 0, align amino acids with the same engine, back-project to a
codon grid. Result (same datasets): SPS ≈ 0.95 / 0.93 at pfs = 0 / 0.2,
phase purity 1.0 — competitive with MAFFT/MACSE at a fraction of runtime.
`refine = TRUE` can apply codon_refine on top but is OFF by default
(it degrades clean high-divergence datasets: 0.954 → 0.871 observed).
Tools recorded separately in the TSV: `msa_codon_protein`,
`msa_codon_protein_refine`.
