# Codon-aware MSA: hare cytochrome b

The `lepus-cytb` vignette used bounded edit distance to check which
records really are cytochrome b. The next step in any marker pipeline is
the multiple sequence alignment — and for a *coding* gene the choice of
aligner mode is a biological decision, not a cosmetic one:

> Frameshift-indifferent alignment inserts indels anywhere, breaking the
> reading frame and turning a clean indel history into apparent
> nonsense. Codon-aware alignment only ever inserts whole codons, so a
> real in-frame deletion looks like what it is.

[`msa_codon()`](https://alrobles.github.io/genoaligner-r/reference/msa_codon.md)
does the second thing. This vignette aligns the 37 real cytb records
from
[`?lepus_cytb`](https://alrobles.github.io/genoaligner-r/reference/lepus_cytb.md)
under the vertebrate mitochondrial code (`gc = 2` — TGA is tryptophan
here, not a stop).

## The alignment

``` r

data(lepus_cytb)

a <- msa_codon(lepus_cytb$sequence, names = lepus_cytb$accession, gc = 2)
a
#> genoaligner MSA: 37 sequences x 1500 columns (mode=codon, gc=2, refine=FALSE)
#>   qc: 23/37 sequences without stop codons; 14 with partial/frameshift blocks
```

One call, ~2.5 seconds, a 37 × 1500 alignment. Because every indel is a
whole codon, every row stays a multiple of 3:

``` r

stopifnot(ncol(a$aligned) %% 3 == 0)
```

## What the QC already knows

`msa_codon` tokenizes each sequence to codons before aligning, and
reports what it found — the chosen reading `frame`, the `n_stops` count
in that frame, and `partial` (a trailing partial codon or held-out
frameshift blocks). The QC travels with the result; no separate pass is
needed:

``` r

qc <- cbind(lepus_cytb[c("accession", "species")], a$qc)
qc[order(-qc$n_stops, -qc$partial), ][1:8, ]
#>         accession     species frame n_stops partial
#> 30 XM_017337704.3   cuniculus     1      22       0
#> 31 XM_051846530.2   cuniculus     0      13       0
#> 1      AF157465.1 granatensis     2       1       0
#> 2      HQ596476.1 granatensis     0       1       0
#> 4      EU285251.1 granatensis     0       1       0
#> 5      EU285250.1 granatensis     0       1       0
#> 19     LC132661.1     timidus     0       1       0
#> 20     LC132659.1     timidus     0       1       0
```

Perfect cytb records score `n_stops = 0, partial = 0`. The records that
do not are exactly the ones the pairwise vignette flags as suspects —
records whose best frame still contains a stop codon are typically the
mislabelled or non-cytb accessions, and `partial > 0` marks fragments
that do not end on a codon boundary. That is the same QC story MACSE
reports, emitted here as an ordinary data.frame column.

## Frame integrity: codon mode vs plain DNA mode

Align the same records as plain DNA:

``` r

d <- msa_align(lepus_cytb$sequence, names = lepus_cytb$accession)
d
#> genoaligner MSA: 37 sequences x 1399 columns (mode=dna)
```

Both alignments are valid outputs of their own rules — but they answer
different questions. In DNA mode an indel can be any width; in codon
mode it cannot. To see it concretely, take two sequences where one
carries a two-codon deletion:

``` r

# a complete-cds record vs a fragment missing two codons
full <- lepus_cytb$sequence[2]
frag <- paste0(substr(full, 1, 60), substr(full, 67, nchar(full)))

cod <- as.matrix(msa_codon(c(full = full, frag = frag), gc = 2))
dna <- as.matrix(msa_align(c(full = full, frag = frag)))

# first gap position in each alignment
gap <- function(row) min(which(row == "-"))
cat("codon-mode gap starts at column:", gap(cod["frag", ]), "\n")
#> codon-mode gap starts at column: 61
cat("  (multiple of 3 ->", (gap(cod["frag", ]) - 1) %% 3 == 0, ")\n")
#>   (multiple of 3 -> TRUE )
cat("dna-mode   gap starts at column:", gap(dna["frag", ]), "\n")
#> dna-mode   gap starts at column: 61
```

In codon mode the deletion is placed as codons — the aligned fragment
reads as clean amino acids end to end. That is the difference that
matters when the alignment feeds a partition file, a codon model, or a
tree search.

## Frameshift refinement

`msa_codon(..., refine = TRUE)` adds a MACSE-style stage-2 pass that
reintroduces real 1-2 nt frameshift events as narrow columns instead of
forcing them through whole-codon indels. The result is a nucleotide MSA
that is *not* guaranteed to be a multiple of 3 — on purpose, because
sequences that genuinely carry frameshifts should not be aligned as if
they did not:

``` r

r <- msa_codon(c("ATGATAATCACCTGA", "ATGATTACACCTGA"), gc = 2, refine = TRUE)
as.character(r)
#> [1] ">seq1\nATG--ATAATCACCTGA" ">seq2\n---ATGATTACACCTGA"
```

## Interoperability

The result is an ordinary character matrix underneath:

``` r

m <- as.matrix(a)
m[1:2, 1:24]
#>            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> AF157465.1 "-"  "-"  "-"  "-"  "-"  "-"  "-"  "-"  "-"  "-"   "-"   "-"   "-"  
#> HQ596476.1 "A"  "T"  "G"  "A"  "C"  "C"  "A"  "A"  "C"  "A"   "T"   "T"   "C"  
#>            [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#> AF157465.1 "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"   "-"  
#> HQ596476.1 "G"   "T"   "A"   "A"   "A"   "A"   "C"   "A"   "C"   "A"   "C"
```

`as.character(a)` emits FASTA records ready for
[`writeLines()`](https://rdrr.io/r/base/writeLines.html), and with
**ape** installed `as_dnabin(a)` coerces straight to `DNAbin`:

``` r

library(ape)
bin <- as_dnabin(a)
bin[1:3, 1:15]
#> 3 DNA sequences in binary format stored in a matrix.
#> 
#> All sequences of same length: 15 
#> 
#> Labels:
#> AF157465.1
#> HQ596476.1
#> MZ203096.1
#> 
#> Base composition:
#>     a     c     g     t 
#> 0.333 0.267 0.133 0.267 
#> (Total: 45 bases)
```

## Where this engine lives

The R package ships the host path of the same C++17 core the GPU library
uses. On the sibling repository the identical engine runs as a
level-batched GPU driver verified bit-exact against this host path — so
an alignment produced here and one produced on the GPU cluster are
byte-identical, not merely similar.
