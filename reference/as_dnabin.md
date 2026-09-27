# Coerce a genoaligner_msa alignment to ape::DNAbin

Requires the ape package (in `Suggests`; used here only).

## Usage

``` r
as_dnabin(x)
```

## Arguments

- x:

  A `genoaligner_msa` object.

## Value

An [`ape::DNAbin`](https://rdrr.io/pkg/ape/man/DNAbin.html) matrix.

## Examples

``` r
if (requireNamespace("ape", quietly = TRUE)) {
    as_dnabin(msa_align(c("ACGTACGT", "ACGTTCGT")))
}
#> 2 DNA sequences in binary format stored in a matrix.
#> 
#> All sequences of same length: 8 
#> 
#> Labels:
#> seq1
#> seq2
#> 
#> Base composition:
#>     a     c     g     t 
#> 0.188 0.250 0.250 0.312 
#> (Total: 16 bases)
```
