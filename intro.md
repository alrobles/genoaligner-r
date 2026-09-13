# NA

`genoaligner` aligns pairs of sequences. It is a **pairwise** aligner:
it compares two sequences you already have, it does not build a multiple
alignment and it does not recruit candidates.

Two entry points cover the two standard pairwise problems:

- [`align()`](https://alrobles.github.io/genoaligner-r/reference/align.md)
  — **edit distance** (Levenshtein; the WFA-equivalent score) with an
  optional CIGAR.
- [`align_sw()`](https://alrobles.github.io/genoaligner-r/reference/align_sw.md)
  — **Smith-Waterman** local alignment with affine gaps, CIGAR, and the
  aligned span’s coordinates.

## Edit distance

``` r

align("ACGTACGT", "ACGTTCGT", smax = 8)
#>   score resolved    cigar rescore_ok wellformed_ok
#> 1     1     TRUE MMMMXMMM       TRUE          TRUE
```

By default `smax = 64`: a pair is *resolved* only when its true
equal-cost edit distance is at most `smax`. If it exceeds the bound the
pair comes back `resolved = FALSE` with an `NA` score rather than a
wrong number — size `smax` to your data to avoid silently under-serving.

Both functions are vectorised over the pair columns, so they map
straight into a pipeline:

``` r

pairs <- data.frame(
  query     = c("ACGTACGT", "GGGG",    "TTTTCCCC", "AAAACCCC"),
  reference = c("ACGTTCGT", "ACGT",    "TTTTCGCG", "AAAATTTT")
)
pairs$result <- align(pairs$query, pairs$reference, smax = 8)
pairs
#>      query reference result.score result.resolved result.cigar
#> 1 ACGTACGT  ACGTTCGT            1            TRUE     MMMMXMMM
#> 2     GGGG      ACGT            3            TRUE         XXMX
#> 3 TTTTCCCC  TTTTCGCG            2            TRUE     MMMMMXMX
#> 4 AAAACCCC  AAAATTTT            4            TRUE     MMMMXXXX
#>   result.rescore_ok result.wellformed_ok
#> 1              TRUE                 TRUE
#> 2              TRUE                 TRUE
#> 3              TRUE                 TRUE
#> 4              TRUE                 TRUE
```

## Smith-Waterman local alignment

``` r

align_sw("TTTACGTGTT", "ACGTGT", scoring = c(2, -3, 5, 2))
#>   score start_i start_j end_i end_j  cigar rescore_ok wellformed_ok too_large
#> 1    12       3       0     8     5 MMMMMM       TRUE          TRUE     FALSE
```

`scoring` is `(match, mismatch, gap_open, gap_extend)`; penalties are
positive and subtracted. The result includes the 0-based coordinates of
the aligned span and the CIGAR over that span.

## The CIGAR convention

Identical for both entry points: `M` and `X` consume a text and a
pattern base (`M` on equal, `X` on differing characters); `I` consumes
an extra **text** base; `D` consumes an extra **pattern** base. Empty
inputs and `score == 0` (no positive-scoring local alignment) are
reported with `NA`/`-1` coordinates.

## GPU backend

The workhorse of this package is a self-contained CPU core that needs no
GPU toolchain. The wavefront/GPU backend (ROCm and CUDA) computes the
same results with acceleration and is tracked in the sibling C++
`genoaligner` repository.
