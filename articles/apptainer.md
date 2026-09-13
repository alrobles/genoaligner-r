# Reproducing the build with Apptainer

The CRAN/CPU core of `genoaligner` builds and runs anywhere R does. When
you want a **reproducible, self-contained environment** – pinned R, Rcpp
and `genoaligner`, decoupled from whatever happens to be installed on
the host – you can build it as an [Apptainer](https://apptainer.org)
(formerly Singularity) container with the recipe that ships in this
repository.

The recipe lives at `containers/genoaligner-r-compute.def`. it is kept
in the repo (and out of the CRAN tarball) so the exact construction
recipe is version-controlled next to the code it pins.

## What you get

A single `.sif` file holding R 4.4.1, Rcpp and `genoaligner`, plus the
build recipe – so a colleague (or a cluster) can reproduce the
environment without trusting your machine’s state.

## The recipe

The file `containers/genoaligner-r-compute.def`:

``` sh
Bootstrap: docker
From: rocker/r-ver:4.4.1

%files
    genoaligner_1.0.0.tar.gz /opt/

%post
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends g++ make ca-certificates
    R -e 'install.packages("Rcpp", repos="https://cloud.r-project.org", quiet=TRUE)'
    R CMD INSTALL /opt/genoaligner_1.0.0.tar.gz
    R -e 'library(genoaligner); stopifnot(align("ACGT","ACGT",smax=8)$score==0); cat("BUILD_CORE_OK\n")'

%runscript
    exec Rscript "$@"
```
