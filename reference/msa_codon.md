# Codon-aware multiple sequence alignment

Codon-aware MSA of coding sequences: each sequence is tokenized to
codons (NCBI genetic-code table `gc`) and aligned over a codon
substitution matrix (BLOSUM62 on the translated amino acids plus a
per-identical-nucleotide bonus and a stop-codon penalty). With the
default `refine = FALSE` every indel is a whole codon, so every aligned
row stays a multiple of 3 and the reading frame is preserved by
construction.

## Usage

``` r
msa_codon(
  seqs,
  names = NULL,
  gc = 1L,
  refine = FALSE,
  iter_refine = 0L,
  fft_band = 0L,
  fft_lags = 4L,
  fft_min_rel = 0.1
)
```

## Arguments

- seqs:

  Character vector; unaligned sequences (no `"-"`). IUPAC nucleotide
  letters in `dna` mode; amino-acid letters in `protein` mode.

- names:

  Optional character vector of sequence names, same length as `seqs`;
  defaults to `names(seqs)` or generated ids. Used for row names of the
  aligned matrix and FASTA output.

- gc:

  Integer; NCBI genetic code: `1L` (standard, default) or `2L`
  (vertebrate mitochondrial — the right table for CYTB/COI/ND).

- refine:

  Logical; run the stage-2 frameshift refinement pass.

- iter_refine:

  Integer \\\ge 0\\; iterative-refinement rounds over the guide-tree
  bipartitions (MAFFT `FFT-NS-i` class): every tree edge splits the
  alignment into two induced sub-profiles, the pair is realigned, and
  the candidate is kept only if the sum-of-pairs score strictly
  improves. `0` (default) keeps the plain progressive alignment. Two
  rounds is a good first value; the pass converges and a further round
  reports no accepted changes.

- fft_band:

  Integer \\\ge 0\\; FFT anchor band half-width in columns. When \\\>
  0\\ (and `iter_refine > 0`) the engine detects homology anchors by
  Fourier cross-correlation of the two profile property signals and runs
  the realignment DP inside a band around them — the literal FFT of
  MAFFT `FFT-NS-i`. Weak anchor evidence falls back to the full DP, and
  the strict-improvement gate applies either way, so this is a
  speed/robustness knob, not a quality risk. `0` (default) uses the full
  matrix.

- fft_lags:

  Integer \\\ge 1\\; number of anchor diagonals whose bands are unioned
  (default `4`). Only used when `fft_band > 0`.

- fft_min_rel:

  Numeric \\\ge 0\\; minimum cosine similarity of the best anchor before
  the band is trusted (default `0.10`). Higher values make the engine
  fall back to the full DP more often.

## Value

A `genoaligner_msa` object. In codon mode `qc` carries measured values
per input sequence: `frame` (chosen reading frame, 0-2), `n_stops` (stop
codons in the tokenized stream) and `partial` (partial/ambiguous codons
and held-out frameshift blocks). A clean coding sequence reports
`n_stops = 0`.

## Details

With `refine = TRUE` a MACSE-style stage-2 pass reintroduces internal
1-2 nt frameshift events as narrow columns, so the output is a
nucleotide MSA no longer guaranteed to be a multiple of 3 — that is the
intended semantics for sequences that truly carry frameshifts.

## Examples

``` r
msa_codon(c("ATGATAATCACC", "ATGATTATCACCTGA"))
#> genoaligner MSA: 2 sequences x 21 columns (mode=codon, gc=1, refine=FALSE)
#>   qc: 2/2 sequences without stop codons; 1 with partial/frameshift blocks
msa_codon(c("ATGATAATCACC", "ATGATTATCACCTGA"), gc = 2)
#> genoaligner MSA: 2 sequences x 15 columns (mode=codon, gc=2, refine=FALSE)
#>   qc: 2/2 sequences without stop codons
```
