# Align pairs by edit distance

Computes the Levenshtein edit distance (substitution = 1, insertion = 1,
deletion = 1) between each pair of sequences, with an optional CIGAR
reconstruction. This is the same distance the WFA formulation computes;
wavefront/GPU speed is a backend concern, the result is identical.

## Usage

``` r
align(pattern, text, smax = 64, with_cigar = TRUE)
```

## Arguments

- pattern:

  Character vector; query sequences (recycled to match `text`).

- text:

  Character vector; reference sequences (recycled to match `pattern`).
  `NA` entries produce `NA` rows.

- smax:

  Numeric; maximum edit distance searched, in `[0, 511]`. Pairs whose
  true distance exceeds it are unresolved by design.

- with_cigar:

  Logical; if `FALSE` only scores are computed (cheaper, no traceback)
  and `cigar` is `NA`.

## Value

A data.frame with one row per pair and columns: `score` (edit distance,
`NA` if unresolved), `resolved` (logical), `cigar` (over {M, X, I, D},
`NA` if unresolved or `!with_cigar`), `rescore_ok` and `wellformed_ok`
(validation flags computed in the library).

CIGAR convention (identical for `align` and `align_sw`): `M`/`X` consume
a text and a pattern base; `I` consumes an extra *text* base; `D`
consumes an extra *pattern* base.

## Details

A pair is *resolved* only when its true distance is at most `smax`;
otherwise `score` is `NA` and the pair is reported as
`resolved = FALSE`. The count of resolved pairs is deliberately visible:
a caller that only reads scores can silently miss that most input was
abandoned when `smax` is set too tightly.

## Examples

``` r
align("ACGTACGT", "ACGTTCGT", smax = 8)
#>   score resolved    cigar rescore_ok wellformed_ok
#> 1     1     TRUE MMMMXMMM       TRUE          TRUE
align(c("AAAACCC", "ACGT"), c("AAAATCC", "ACCT"), smax = 8)
#>   score resolved   cigar rescore_ok wellformed_ok
#> 1     1     TRUE MMMMXMM       TRUE          TRUE
#> 2     1     TRUE    MMXM       TRUE          TRUE
```
