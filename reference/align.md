# Align pairs by edit distance (deprecated alias)

`align()` is deprecated in favour of
[`align_edit()`](https://alrobles.github.io/genoaligner-r/reference/align_edit.md):
the generic name collides with Biostrings/IRanges exports in sessions
that load both. It still works and returns identical results; new code
should call
[`align_edit()`](https://alrobles.github.io/genoaligner-r/reference/align_edit.md).

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

See
[`align_edit`](https://alrobles.github.io/genoaligner-r/reference/align_edit.md).
