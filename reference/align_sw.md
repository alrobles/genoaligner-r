# Smith-Waterman local alignment

Local alignment with affine gap penalties {match, mismatch, gap_open,
gap_extend}. Penalties are positive and subtracted. A pair with
`score == 0` has no positive-scoring local alignment (disjoint
sequences); all four coordinates are `-1` and `cigar` is empty.

## Usage

``` r
align_sw(text, pattern, scoring = c(2L, -3L, 5L, 2L), with_cigar = TRUE)
```

## Arguments

- text:

  Character vector; row-axis sequences (recycled).

- pattern:

  Character vector; column-axis sequences (recycled).

- scoring:

  Numeric vector of length 4, named or positional:
  `(match, mismatch, gap_open, gap_extend)`. Recycled if a single vector
  is given. The defaults are an example, not a recommendation.

- with_cigar:

  Logical; if `FALSE`, only score and end coordinates are computed
  (start coordinates and `cigar` are `NA`/`-1`).

## Value

A data.frame with one row per pair and columns: `score`, `start_i`,
`start_j`, `end_i`, `end_j` (0-based coordinates of the aligned span,
`-1` when score is 0/unknown), `cigar` (over the span), `rescore_ok`,
`wellformed_ok` and `too_large` (the pair's matrix exceeded the
supported size limit).

## Examples

``` r
align_sw("TTTACGTGTT", "ACGTGT", scoring = c(2, -3, 5, 2))
#>   score start_i start_j end_i end_j  cigar rescore_ok wellformed_ok too_large
#> 1    12       3       0     8     5 MMMMMM       TRUE          TRUE     FALSE
align_sw(c("ACGTACGT", "GGGG"), c("ACGTTCGT", "ACGT"), c(2, -3, 5, 2))
#>   score start_i start_j end_i end_j    cigar rescore_ok wellformed_ok too_large
#> 1    11       0       0     7     7 MMMMXMMM       TRUE          TRUE     FALSE
#> 2     2       0       2     0     2        M       TRUE          TRUE     FALSE
```
