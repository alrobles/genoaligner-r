## R4 diagnostics: guided re-tokenization evidence used by
## refine="auto". guide_bonus_cpp maps each raw nucleotide to its
## pass-1 alignment column and returns a leave-one-out agreement bonus;
## codon_encode_blocks_cpp re-tiles each sequence (codon blocks of 3,
## held-out blocks of -d) and reports block provenance: internal
## frameshift blocks (codons on both sides) vs terminal raggedness
## (lead_nt / trail_nt). The replication is parity-checked against the
## vendored genomsa::codon_encode_guided at default parameters.

.aln_rows <- function(a) apply(as.matrix(a), 1L, paste0, collapse = "")

test_that("guide_bonus: per-row bonus vector has one entry per nucleotide", {
    sq <- c(a = "ATGAAACCCGGG", b = "ATGAAACCCGGGAAA",
            c = "ATGCCCGGGTAA")
    rows <- .aln_rows(msa_codon(unname(sq), backbone = "protein"))
    g <- genoaligner:::guide_bonus_cpp(unname(sq), rows)
    expect_length(g, 3L)
    expect_true(all(vapply(g, length, 0L) == nchar(unname(sq))))
    expect_true(all(vapply(g, function(v) all(is.finite(v)), TRUE)))
})

test_that("encode_blocks: nucleotide accounting is exact", {
    sq <- c(a = "ATGAAACCCGGGTTT",      # clean 15
            b = "ATGAAAGCCCGGGTTT",    # +1 nt internal insertion
            c = "ATGAAACCCGGGTT")      # trailing partial (14)
    eb <- genoaligner:::codon_encode_blocks_cpp(
        unname(sq), 1L, 0, replicate(3L, rep(0, 1), simplify = FALSE), -1)
    accounted <- eb$lead_nt + 3L * eb$n_codons + eb$fs_nt + eb$trail_nt
    expect_equal(accounted, nchar(unname(sq)))
    # the ragged rows must hold out their extra nt in *terminal* blocks
    expect_equal(eb$lead_nt[3L] + eb$trail_nt[3L], 2L)
    expect_gte(eb$lead_nt[2L] + eb$trail_nt[2L], 1L)
})

test_that("encode_blocks: parity with vendored codon_encode_guided", {
    sq <- c(a = "ATGAAACCCGGGTTT", b = "ATGAAAGCCCGGGTTT",
            c = "ATGAAACCCGGGTT", d = "ATGCCCGGGTAAATA")
    rows <- .aln_rows(msa_codon(unname(sq), backbone = "protein"))
    g <- genoaligner:::guide_bonus_cpp(unname(sq), rows)
    eb <- genoaligner:::codon_encode_blocks_cpp(unname(sq), 1L, 40, g, -1)
    expect_true(all(eb$parity == 1L))
    expect_true(all(eb$vend_partial >= 0L))
    # internal blocks are a subset of all held-out blocks
    expect_true(all(eb$n_fs <= eb$vend_partial))
})

test_that("encode_blocks: fs_enc override suppresses internal blocks", {
    # identical rows cannot tile a sequence needing internal evidence;
    # a strongly supported internal block should survive a raised cost
    sq <- c(a = paste0("ATG", strrep("GCA", 12L), "TAA"),
            b = paste0("ATG", strrep("GCA", 12L), "CTAGCTAA"),
            c = paste0("ATG", strrep("GCA", 12L), "CTAGCTAATAA"))
    rows <- .aln_rows(msa_codon(unname(sq), backbone = "protein"))
    g <- genoaligner:::guide_bonus_cpp(unname(sq), rows)
    lo <- genoaligner:::codon_encode_blocks_cpp(unname(sq), 1L, 40, g, 80)
    hi <- genoaligner:::codon_encode_blocks_cpp(unname(sq), 1L, 40, g, 1e6)
    expect_true(all(hi$n_fs <= lo$n_fs))          # monotone in cost
    accounted <- hi$lead_nt + 3L * hi$n_codons + hi$fs_nt + hi$trail_nt
    expect_equal(accounted, nchar(unname(sq)))    # accounting still exact
})

test_that("refine='auto' is phase-neutral on ragged CDS with no internal
          frameshift evidence", {
    # five rows, three with terminal partial codons, all in-frame content
    sq <- c(a = paste0("ATG", strrep("GCA", 15L)),
            b = paste0("ATG", strrep("GCA", 15L), "G"),
            c = paste0("ATG", strrep("GCA", 15L), "GC"),
            d = paste0("AT", strrep("GCA", 15L)),
            e = paste0("ATG", strrep("GCA", 15L), "GCAT"))
    a0 <- msa_codon(sq, backbone = "protein")
    aA <- msa_codon(sq, backbone = "protein", refine = "auto")
    expect_identical(as.matrix(a0), as.matrix(aA))
})

test_that("selective refine honours the flag mask (plumbing)", {
    sq <- c(a = "ATGAAACCCGGGTTT", b = "ATGAAAGCCCGGGTTT",
            c = "ATGAAACCCGGGTT")
    rows <- .aln_rows(msa_codon(unname(sq), backbone = "protein"))
    none <- genoaligner:::codon_refine_sel_cpp(
        unname(sq), rows, rep(FALSE, 3L), 1L, 1L)
    # no flags: no sequence content changes (trailing empty codon slots
    # may still be trimmed by the merge)
    expect_equal(gsub("-", "", none$aln), unname(sq))
    sel <- genoaligner:::codon_refine_sel_cpp(
        unname(sq), rows, c(FALSE, TRUE, TRUE), 1L, 1L)
    expect_true(all(nchar(sel$aln) == nchar(sel$aln[1])))
    expect_equal(gsub("-", "", sel$aln), unname(sq))
    # the flag helper itself returns one logical per row
    fl <- genoaligner:::auto_fs_flag(sq, rows, 1L)
    expect_identical(length(fl), 3L)
    expect_type(fl, "logical")
})
