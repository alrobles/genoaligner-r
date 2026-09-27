#' Cytochrome b sequences from four leporid species
#'
#' Real mitochondrial cytochrome b (cytb) records downloaded from NCBI
#' nuccore for three hare species (*Lepus granatensis*, *L. timidus*,
#' *L. europaeus*) — all flagged as population-genetic candidates in the
#' EcoSeek genominer database (\url{https://genominer.ecoseek.org}) — plus
#' the European rabbit (*Oryctolagus cuniculus*) as an outgroup.
#'
#' Sequences are shipped exactly as retrieved: mixed lengths, untrimmed,
#' and including two mislabelled records (XM_017337704.3 and XM_051846530.2,
#' predicted nuclear mRNAs that a "cytochrome b" title search returns but
#' that are not the mitochondrial gene). The `lepus-cytb` vignette uses them
#' to show that a near-zero local alignment score flags a mislabelled record.
#'
#' @format A data.frame with 37 rows and 4 columns:
#' \describe{
#'   \item{accession}{NCBI nuccore accession (character).}
#'   \item{species}{Species epithet: \code{granatensis}, \code{timidus},
#'     \code{europaeus} or \code{cuniculus} (character).}
#'   \item{title}{Full NCBI FASTA header (character).}
#'   \item{sequence}{Nucleotide sequence as downloaded, untrimmed
#'     (character, 400--1200 bp; MN098958.1 contains one IUPAC ambiguity
#'     code, W).}
#' }
#'
#' @source NCBI nuccore, downloaded 2026-09-24; see
#'   `data-raw/lepus_cytb.R` in the source repository for the exact
#'   eutils queries and snapshot FASTAs.
#'
#' @examples
#' data(lepus_cytb)
#' table(lepus_cytb$species)
#' align_edit(lepus_cytb$sequence[1:3], lepus_cytb$sequence[2], smax = 200)
"lepus_cytb"
