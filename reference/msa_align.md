# Multiple sequence alignment (DNA or protein)

Progressive MSA: fragment-corrected k-mer distances, a deterministic
neighbor-joining guide tree, post-order profile-vs-profile semiglobal
alignment (free end gaps, position-specific gap penalties). The pipeline
is deterministic: identical input always produces identical output.

## Usage

``` r
msa_align(seqs, names = NULL, mode = c("dna", "protein"), guide = "nj")
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
