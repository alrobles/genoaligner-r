# genoaligner: GPU-portable pairwise sequence alignment

Pairwise sequence alignment from one portable C++17 core: edit distance
(Levenshtein / WFA-equivalent) and Smith-Waterman local alignment, both
with score and CIGAR reconstruction. The core builds and runs anywhere
(it is the CPU backend, so it needs no GPU toolchain); the sibling C++
library tracks the GPU (ROCm/CUDA) backend that computes the same
results with wavefront acceleration.

## Details

The two entry points are
[`align`](https://alrobles.github.io/genoaligner-r/reference/align.md)
(edit distance bounded by `smax`) and
[`align_sw`](https://alrobles.github.io/genoaligner-r/reference/align_sw.md)
(local Smith-Waterman with affine gaps). Both are vectorised over pairs,
so a two-column pipeline (query, reference) maps straight to a
data.frame of scores and CIGARs.

## See also

Useful links:

- <https://github.com/alrobles/genoaligner-r>

- <https://alrobles.github.io/genoaligner-r/>

- Report bugs at <https://github.com/alrobles/genoaligner-r/issues>

## Author

**Maintainer**: Angel Robles-Fernandez <a.l.robles.fernandez@gmail.com>

Authors:

- Angel Robles-Fernandez <a.l.robles.fernandez@gmail.com>
