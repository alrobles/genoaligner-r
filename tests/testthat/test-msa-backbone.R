## protein-backbone path of msa_codon(): translate -> protein MSA ->
## back-project to a codon grid. These tests pin the contract that makes
#  the route useful: whole-codon indels, equal row lengths, frame 0
#  translation, honest QC for stops/partials.

test_that("codon_tr_frame: standard code (gc=1) translation", {
    tab <- genoaligner:::.codon_table$gc1
    tr <- genoaligner:::codon_tr_frame("ATGAAATTTTGA", tab)   # M K F *
    expect_equal(tr$aa, "MKFX")
    expect_equal(tr$n_stops, 1L)
    expect_equal(tr$raw, c("ATG", "AAA", "TTT", "TGA"))
    expect_equal(nchar(tr$tail), 0L)
})

test_that("codon_tr_frame: vertebrate mitochondrial code (gc=2)", {
    tab <- genoaligner:::.codon_table$gc2
    # TGA is Trp (not stop), ATA is Met, AGA is a stop under gc=2
    tr <- genoaligner:::codon_tr_frame("ATGTGAATAAGA", tab)   # M W M *
    expect_equal(tr$aa, "MWMX")
    expect_equal(tr$n_stops, 1L)
    tr1 <- genoaligner:::codon_tr_frame("ATGTGAATAAGA", genoaligner:::.codon_table$gc1)
    expect_equal(tr1$aa, "MXIR")                # gc=1: TGA stops, AGA is Arg
    expect_equal(tr1$n_stops, 1L)
})

test_that("backbone=protein: rectangular output, names preserved", {
    s <- c(a = "ATGAAATTT", b = "ATGAAATTTGGG", c = "ATGCCCTTT")
    a <- msa_codon(s, backbone = "protein")
    m <- as.matrix(a)
    expect_equal(rownames(m), c("a", "b", "c"))
    expect_true(all(apply(m, 1, function(r) nchar(paste0(r, collapse = "")))
                    == ncol(m)))
    expect_equal(a$params$backbone, "protein")
})

test_that("backbone=protein: all in-frame input yields a multiple-of-3 width", {
    s <- c("ATGAAATTT", "ATGAAATTTGGG", "ATGCCCTTT")
    a <- msa_codon(s, backbone = "protein")
    expect_equal(ncol(as.matrix(a)) %% 3L, 0L)
})

test_that("backbone=protein: protein gaps become exactly three nt gaps", {
    s <- c(a = "ATGAAATTT", b = "ATGAAATTTGGG")
    a <- msa_codon(s, backbone = "protein")
    m <- as.matrix(a)
    # seq a is one codon shorter: its gap run(s) must be a multiple of 3
    runs <- rle(m["a", ])
    gaplens <- runs$lengths[runs$values == "-"]
    expect_true(all(gaplens %% 3L == 0L))
})

test_that("backbone=protein: codon phase purity is 1 on in-frame input", {
    s <- c("ATGAAATTTCCC", "ATGAAAGGGCCC", "ATGCCCTTTAAA")
    a <- msa_codon(s, backbone = "protein")
    m <- as.matrix(a)
    ph <- t(vapply(seq_len(nrow(m)), function(i) {
        b <- m[i, ] != "-"; (cumsum(b) - 1L) %% 3L
    }, integer(ncol(m))))
    pure <- vapply(seq_len(ncol(m)), function(j)
        length(unique(ph[m[, j] != "-", j])) <= 1L, logical(1))
    expect_true(all(pure))
})

test_that("backbone=protein: non-mult-3 seq keeps dangling nt as end insertions", {
    s <- c(a = "ATGAAATT", b = "ATGAAATTT")      # a: 8 nt (2-nt tail)
    a <- msa_codon(s, backbone = "protein")
    m <- as.matrix(a)
    expect_equal(a$qc$partial[a$qc$partial > 0], 1)
    # the 2 dangling nt appear verbatim somewhere in row a's tail block
    expect_match(paste0(m["a", ], collapse = ""), "TT-*$")
})

test_that("backbone=protein: ambiguous codons are placed, not dropped", {
    s <- c(a = "ATGNNNTTT", b = "ATGAAATTT")
    a <- msa_codon(s, backbone = "protein")
    m <- as.matrix(a)
    expect_match(paste0(m["a", ], collapse = ""), "NNN")
    expect_gte(a$qc$partial[1], 1L)
})

test_that("backbone=codon default is unchanged", {
    s <- c(a = "ATGAAATTT", b = "ATGAAATTTGGG")
    expect_equal(msa_codon(s)$params$backbone, "codon")
    expect_equal(msa_codon(s, backbone = "codon")$params$backbone, "codon")
})

test_that("invalid backbone is rejected", {
    expect_error(msa_codon("ATGAAA", backbone = "x"))
})

test_that("backbone=protein + refine=TRUE runs the stage-2 pass", {
    s <- c("ATGATAATCACCTGA", "ATGATTACACCTGA")
    a <- msa_codon(s, gc = 2, backbone = "protein", refine = TRUE)
    expect_s3_class(a, "genoaligner_msa")
    expect_true(a$params$refine)
})

test_that("refine=auto: in-frame input is returned exactly as the backbone", {
    s <- c(a = "ATGAAATTTCCC", b = "ATGAAAGGGCCC", c = "ATGCCCTTTAAA")
    auto <- msa_codon(s, backbone = "protein", refine = "auto")
    plain <- msa_codon(s, backbone = "protein")
    expect_identical(as.matrix(auto), as.matrix(plain))
    # same invariance on the codon backbone
    auto_c <- msa_codon(s, refine = "auto")
    plain_c <- msa_codon(s)
    expect_identical(as.matrix(auto_c), as.matrix(plain_c))
})

test_that("refine=auto: a non-mult-3 sequence can emit event columns", {
    s <- c(a = "ATGAAATTTCC", b = "ATGAAATTTCCC")   # a: 11 nt
    a <- msa_codon(s, backbone = "protein", refine = "auto")
    expect_s3_class(a, "genoaligner_msa")
    expect_equal(a$params$refine, "auto")
    ac <- msa_codon(s, refine = "auto")
    expect_s3_class(ac, "genoaligner_msa")
    # every row keeps the same length
    m <- as.matrix(ac)
    expect_true(all(apply(m, 1, function(r) nchar(paste0(r, collapse = "")))
                    == ncol(m)))
})

test_that("refine rejects invalid values", {
    expect_error(msa_codon("ATGAAA", refine = "x"))
    expect_error(msa_codon("ATGAAA", backbone = "protein", refine = "x"))
})
