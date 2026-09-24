# Build the `lepus_cytb` package dataset.
#
# Real cytochrome b (cytb) records downloaded from NCBI nuccore on 2026-09-24
# for three hare species that appear as population-genetic candidates in the
# EcoSeek genominer database (https://genominer.ecoseek.org), plus the European
# rabbit as an outgroup. The snapshot FASTAs live in data-raw/fasta/ so the
# dataset can be rebuilt without network access.
#
# Queries used (db=nuccore, rettype=fasta):
#   "Lepus granatensis"[organism] AND (cytb[title] OR "cytochrome b"[title])
#       AND 400:1200[slen]                          -> 13 records
#   "Lepus timidus"[organism]    AND "cytochrome b"[title]
#       AND 1000:1200[slen], sort=relevance          -> 8 records
#   "Lepus europaeus"[organism]  AND "cytochrome b"[title]
#       AND 1000:1200[slen], sort=relevance          -> 8 records
#   "Oryctolagus cuniculus"[organism] AND "cytochrome b"[title]
#       AND 1000:1200[slen], sort=relevance          -> 8 records
#
# Sequences are shipped exactly as downloaded: mixed lengths and untrimmed,
# on purpose -- the lepus-cytb vignette maps the short fragments onto the
# complete gene with align_sw() and crops to the common window.
#
# Two Oryctolagus records (XM_017337704.3, XM_051846530.2) are predicted
# *nuclear* mRNAs that NCBI returns under a "cytochrome b" title search but
# that are not the mitochondrial gene. They are kept in the dataset on
# purpose: the vignette uses them to show that a near-zero local score flags
# a mislabelled record.

read_fasta <- function(path) {
    l <- readLines(path)
    hdr <- grep("^>", l)
    data.frame(
        accession = sub(" .*", "", sub("^>", "", l[hdr])),
        title = sub("^>", "", l[hdr]),
        sequence = vapply(seq_along(hdr), function(i) {
            end <- if (i < length(hdr)) hdr[i + 1L] - 1L else length(l)
            paste0(l[(hdr[i] + 1L):end], collapse = "")
        }, character(1)),
        stringsAsFactors = FALSE
    )
}

files <- c(
    granatensis = "data-raw/fasta/Lepus_granatensis.fasta",
    timidus     = "data-raw/fasta/Lepus_timidus.fasta",
    europaeus   = "data-raw/fasta/Lepus_europaeus.fasta",
    cuniculus   = "data-raw/fasta/Oryctolagus_cuniculus.fasta"
)

lepus_cytb <- do.call(rbind, Map(function(f, sp) {
    d <- read_fasta(f)
    d$species <- sp
    d
}, files, names(files)))

# drop duplicated records (AF157465.1 was returned by two granatensis queries)
lepus_cytb <- lepus_cytb[!duplicated(lepus_cytb$accession), ]
lepus_cytb <- lepus_cytb[, c("accession", "species", "title", "sequence")]
rownames(lepus_cytb) <- NULL

# MN098958.1 carries one IUPAC ambiguity code (W = A|T); kept as downloaded.
stopifnot(
    all(nchar(lepus_cytb$sequence) > 0L),
    all(grepl("^[ACGTNWRYKMSBDHV]+$", lepus_cytb$sequence)),
    nrow(lepus_cytb) == 37L
)

usethis::use_data(lepus_cytb, overwrite = TRUE)
