# Multiple sequence alignment (DNA or protein)

Progressive MSA: fragment-corrected k-mer distances, a deterministic
neighbor-joining guide tree, post-order profile-vs-profile semiglobal
alignment (free end gaps, position-specific gap penalties). The pipeline
is deterministic: identical input always produces identical output.

## Usage

``` r
msa_align(
  seqs,
  names = NULL,
  mode = c("dna", "protein"),
  guide = "nj",
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

- mode:

  `"dna"` (default) or `"protein"`.

- guide:

  Guide-tree construction; only `"nj"` (deterministic neighbor joining)
  exists today.

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

An object of class `genoaligner_msa` with components `aligned`
(character matrix, one row per input sequence, one column per alignment
position), `qc` (per-sequence data.frame; all `NA` outside codon mode)
and `params` (echo of the effective request). See
[`msa_codon`](https://alrobles.github.io/genoaligner-r/reference/msa_codon.md)
for the codon-aware mode.

## Examples

``` r
msa_align(c("ACGTACGTAA", "ACGTTCGT", "ACGTACG"))
#> genoaligner MSA: 3 sequences x 10 columns (mode=dna)
msa_align(c(insulin = "MALWMRLLPLL", insulin2 = "MALWTRLLPLL"),
          mode = "protein")
#> genoaligner MSA: 2 sequences x 11 columns (mode=protein)
```
